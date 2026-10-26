// yruler: Y satellite array lengths from reads
// count: SAM on stdin -> marker k-mer counts per class
// estimate: class sums -> array lengths, Y dosage
use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicU32, Ordering::Relaxed};
use std::sync::{mpsc, Arc, Mutex};
use std::time::Instant;

const K: usize = 31;
const MASK: u64 = (1u64 << 62) - 1;
const EMPTY: u64 = u64::MAX;
const CHUNK: usize = 8 << 20;

static LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'C' as usize] = 1; t[b'G' as usize] = 2; t[b'T' as usize] = 3;
    t[b'a' as usize] = 0; t[b'c' as usize] = 1; t[b'g' as usize] = 2; t[b't' as usize] = 3;
    t
};


#[inline(always)]
fn hk(x: u64) -> u64 { let m = x.wrapping_mul(0x9e3779b97f4a7c15); m ^ (m >> 29) }   // k-mer hash: one multiply, high bits folded down

fn die(m: &str) -> ! { eprintln!("ERROR yruler: {m}"); std::process::exit(1) }

// word Bloom filter, checked only on k-mer cache misses: 4 bits in one 64-bit word per key (one load, one AND)
struct Bloom { b: Vec<u64>, mask: u64 }
impl Bloom {
    fn new(n: usize, bits_per_key: usize) -> Self {
        let nw = ((n * bits_per_key + 63) / 64).next_power_of_two();
        let mut b: Vec<u64> = Vec::with_capacity(nw); hugepage(&b); b.resize(nw, 0); Bloom { b, mask: nw as u64 - 1 }
    }
    #[inline(always)]
    fn pat(h: u64) -> u64 { (1u64 << (h & 63)) | (1u64 << ((h >> 6) & 63)) | (1u64 << ((h >> 12) & 63)) | (1u64 << ((h >> 18) & 63)) }
    #[inline(always)]
    fn at(&self, h: u64) -> *const u64 { unsafe { self.b.as_ptr().add(((h >> 32) & self.mask) as usize) } }
    fn insert(&mut self, h: u64) { let i = ((h >> 32) & self.mask) as usize; self.b[i] |= Self::pat(h); }
    #[inline(always)]
    fn has(&self, h: u64) -> bool { let m = Self::pat(h); unsafe { *self.at(h) & m == m } }
}

// open-addressing table, one 16-byte slot per k-mer: key, count, marker index | class bits << 25 (one cache line per hit)
#[repr(C, align(16))]
struct Slot { key: u64, cnt: AtomicU32, meta: u32 }
const IDX: u32 = (1 << 25) - 1;
fn hugepage<T>(v: &Vec<T>) {   // ask the kernel for transparent huge pages before the first touch (fewer TLB misses)
    let (p, n) = (v.as_ptr() as usize, v.capacity() * std::mem::size_of::<T>()); let a = (p + (1 << 21) - 1) & !((1 << 21) - 1);
    if n > (4 << 21) { unsafe { libc::madvise(a as *mut libc::c_void, (p + n - a) & !((1 << 21) - 1), libc::MADV_HUGEPAGE); } }
}
struct Table { s: Vec<Slot>, mask: usize }
impl Table {
    fn new(n: usize, thr: usize) -> Self {   // slots written in parallel (first touch spread over threads)
        let s = (n * 2).next_power_of_two(); let mut v: Vec<Slot> = Vec::with_capacity(s); hugepage(&v);
        let (p, per) = (v.as_mut_ptr() as usize, s.div_ceil(thr));
        std::thread::scope(|sc| for t in 0..thr { sc.spawn(move || { let p = p as *mut Slot; for i in t * per..((t + 1) * per).min(s) { unsafe { p.add(i).write(Slot { key: EMPTY, cnt: AtomicU32::new(0), meta: 0 }); } } }); });
        unsafe { v.set_len(s); } Table { s: v, mask: s - 1 }
    }
    // returns the marker index of key, inserting it with index val; class bit c is ORed in
    fn insert(&mut self, key: u64, h: u64, val: u32, c: usize) -> (u32, bool) {
        let mut i = (h >> 7) as usize & self.mask;
        loop {
            let sl = &mut self.s[i];
            if sl.key == EMPTY { sl.key = key; sl.meta = val | (1 << (25 + c)); return (val, true); }
            if sl.key == key { sl.meta |= 1 << (25 + c); return (sl.meta & IDX, false); }
            i = (i + 1) & self.mask;
        }
    }
    #[inline(always)]
    fn at(&self, h: u64) -> *const Slot { unsafe { self.s.as_ptr().add((h >> 7) as usize & self.mask) } }
    #[inline(always)]
    fn find(&self, key: u64, h: u64) -> Option<usize> {
        let mut i = (h >> 7) as usize & self.mask;
        loop {
            let sl = &self.s[i];
            if sl.key == key { return Some(i); }
            if sl.key == EMPTY { return None; }
            i = (i + 1) & self.mask;
        }
    }
}
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")] unsafe { std::arch::x86_64::_mm_prefetch(p as *const i8, std::arch::x86_64::_MM_HINT_T0); }
}

fn canon(s: &[u8]) -> Option<u64> {
    if s.len() != K { return None; }
    let (mut f, mut r) = (0u64, 0u64);
    for &b in s { let c = LUT[b as usize]; if c > 3 { return None; } f = (f << 2) | c as u64; r = (r >> 2) | ((3 - c as u64) << 60); }
    Some(f.min(r))
}
fn decode(x: u64) -> [u8; K] { let mut o = [0u8; K]; for i in 0..K { o[i] = b"ACGT"[((x >> (2 * (K - 1 - i))) & 3) as usize]; } o }

struct Markers { key: Vec<u64>, order: Vec<u32>, ocls: Vec<usize>, gc: Vec<u8>, names: Vec<String>, tab: Table }
// class name of a header line: text before the first '_' (trailing whitespace trimmed first)
fn hname(h: &[u8]) -> &[u8] { let mut e = h.len(); while e > 0 && h[e - 1].is_ascii_whitespace() { e -= 1; } let h = &h[..e]; &h[..h.iter().position(|&b| b == b'_').unwrap_or(h.len())] }
// file cut at line ends into thr pieces parsed in parallel, then keys inserted in file order (first occurrence = index)
// masked superstring: uppercase base = a k-mer starts there
fn expand_ms(raw: &[u8]) -> Vec<u8> {
    let mut o = Vec::with_capacity(raw.len() * 33); let mut seq: Vec<u8> = vec![]; let mut hd: &[u8] = b"";
    let mut flush = |hd: &[u8], seq: &mut Vec<u8>, o: &mut Vec<u8>| { if hd.is_empty() { if !seq.is_empty() { die("sequence before the first header"); } return; }
        o.push(b'>'); o.extend_from_slice(hname(hd)); o.push(b'\n');
        for i in 0..(seq.len() + 1).saturating_sub(K) { if seq[i].is_ascii_uppercase() { o.extend(seq[i..i + K].iter().map(|b| b.to_ascii_uppercase())); o.push(b'\n'); } }
        seq.clear(); };
    for l in raw.split(|&b| b == b'\n') { let l = l.strip_suffix(b"\r").unwrap_or(l); if l.first() == Some(&b'>') { flush(hd, &mut seq, &mut o); hd = &l[1..]; } else { seq.extend_from_slice(l); } }
    flush(hd, &mut seq, &mut o); o
}
fn load(path: &str, thr: usize) -> Markers {
    let raw = std::fs::read(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    // one k-mer per line, else masked superstring
    let ms = raw.split(|&b| b == b'\n').find(|l| !l.is_empty() && l[0] != b'>').map_or(false, |l| { let l = l.strip_suffix(b"\r").unwrap_or(l); l.len() != K || l.iter().any(|b| b.is_ascii_lowercase()) });
    let raw = if ms { let t = Instant::now(); let e = expand_ms(&raw); eprintln!("markers: masked superstring {} bytes expanded to {} bytes in {:.1} s", raw.len(), e.len(), t.elapsed().as_secs_f64()); e } else { raw };
    let n = raw.len();
    let mut cuts = vec![0usize]; for t in 1..thr { let mut p = (n * t / thr).max(*cuts.last().unwrap()); while p < n && p > 0 && raw[p - 1] != b'\n' { p += 1; } cuts.push(p); } cuts.push(n);
    // per piece: (key, hash, G+C, local header index or 255 = before the first header of the piece), header names in order
    let parts: Vec<(Vec<(u64, u64, u8, u8)>, Vec<&[u8]>)> = std::thread::scope(|sc| {
        let hs: Vec<_> = cuts.windows(2).map(|w| { let raw = &raw[w[0]..w[1]]; sc.spawn(move || {
            let (mut v, mut nm, mut cur): (Vec<(u64, u64, u8, u8)>, Vec<&[u8]>, u8) = (Vec::with_capacity(raw.len() / 64 + 16), vec![], 255);
            for l in raw.split(|&b| b == b'\n') {
                if l.is_empty() { continue; }
                if l[0] == b'>' { let c = hname(&l[1..]); if nm.last() != Some(&c) { if nm.len() == 254 { die("more than 254 header names in one piece"); } nm.push(c); } cur = (nm.len() - 1) as u8; continue; }
                let l = l.strip_suffix(b"\r").unwrap_or(l);
                let x = canon(l).unwrap_or_else(|| die(&format!("marker line of length {} or non-ACGT", l.len())));
                v.push((x, hk(x), l.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count() as u8, cur));
            }
            (v, nm) }) }).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut names: Vec<String> = vec![]; let mut cur = usize::MAX; let tot: usize = parts.iter().map(|p| p.0.len()).sum();
    if tot > IDX as usize { die("more than 33,554,431 markers"); }
    let mut tab = Table::new(tot, thr); let (mut key, mut order, mut ocls, mut gc) = (vec![], Vec::with_capacity(tot), Vec::with_capacity(tot), Vec::with_capacity(tot));
    for (v, nm) in &parts {
        let gl: Vec<usize> = nm.iter().map(|c| { let c = String::from_utf8_lossy(c).to_string(); names.iter().position(|x| *x == c).unwrap_or_else(|| { names.push(c); names.len() - 1 }) }).collect();
        if names.len() > 7 { die("more than 7 classes (class bits 25-31 of the slot)"); }
        for (k, &(x, h, g, lc)) in v.iter().enumerate() {
            if lc != 255 { cur = gl[lc as usize]; } else if cur == usize::MAX { die("sequence before the first header"); }
            if let Some(nx) = v.get(k + 16) { prefetch(tab.at(nx.1)); }
            let (i, new) = tab.insert(x, h, key.len() as u32, cur); if new { key.push(x); }
            order.push(i); ocls.push(cur); gc.push(g);
        }
        if let Some(&c) = gl.last() { cur = c; }
    }
    Markers { key, order, ocls, gc, names, tab }
}

// per-thread bin record: reads, reads MAPQ >= 20, duplicates, then per class hits, hits MAPQ >= 20, reads with a hit
struct Local { contigs: HashMap<Vec<u8>, u32>, cname: Vec<Vec<u8>>, bins: HashMap<(u32, u64), usize>, rec: Vec<(u32, u64, Vec<u64>)>, reads: u64, kmers: u64, hits: u64, cmiss: u64, bpos: u64, cache: Vec<Ent> }
// per-thread cache of recent k-mers (sorted input: a k-mer comes back in the next reads); counts added to the table on eviction
#[derive(Clone, Copy)]
struct Ent { key: u64, sc: u32, pend: u32 }   // 16 bytes; sc = table slot | class bits << 25, 0 = not a marker (markers have >= 1 class bit)
fn evict(e: &Ent, tab: &Table) { if e.sc != 0 && e.pend > 0 { tab.s[(e.sc & IDX) as usize].cnt.fetch_add(e.pend, Relaxed); } }

#[inline(always)]
fn mc(s: &[u8], c: u8) -> Option<usize> { let p = unsafe { libc::memchr(s.as_ptr() as *const libc::c_void, c as i32, s.len()) }; if p.is_null() { None } else { Some(p as usize - s.as_ptr() as usize) } }   // glibc SIMD memchr
#[inline(always)]
fn field<'a>(l: &'a [u8], q: &mut usize) -> &'a [u8] {   // next tab-separated field from q (last field: to line end)
    if *q > l.len() { die("SAM line with too few fields"); }
    let r = &l[*q..]; let e = mc(r, b'\t').unwrap_or(r.len()); *q += e + 1; &r[..e]
}
fn num(s: &[u8]) -> u64 { let mut v = 0u64; for &b in s { if !b.is_ascii_digit() { die("non-numeric SAM field"); } v = v * 10 + (b - b'0') as u64; } v }

fn work(chunk: &[u8], loc: &mut Local, bloom: &Bloom, tab: &Table, ncl: usize, bin: u64, fa: &mut Vec<u8>, kb: &mut Vec<(u64, u64, bool)>) {
    let w = 3 + 3 * ncl; let mut last: (Vec<u8>, u64, usize) = (Vec::new(), u64::MAX, 0);
    let mut h = vec![0u64; ncl];
    let mut p0 = 0;
    while p0 < chunk.len() {
        let e0 = mc(&chunk[p0..], b'\n').map_or(chunk.len(), |i| p0 + i); let line = &chunk[p0..e0]; p0 = e0 + 1;
        if line.is_empty() || line[0] == b'@' { continue; }
        let mut q = 0; let it = &mut q;
        field(line, it); let flag = field(line, it); let rname = field(line, it); let pos = num(field(line, it)); let mapq = field(line, it);
        for _ in 0..4 { field(line, it); } let seq = field(line, it);
        let mq = num(mapq) >= 20; let fl = num(flag); let b = pos / bin;
        let ri = if last.1 == b && last.0 == rname { last.2 } else {
            let ci = match loc.contigs.get(rname) { Some(&c) => c, None => { let c = loc.cname.len() as u32; loc.contigs.insert(rname.to_vec(), c); loc.cname.push(rname.to_vec()); c } };
            let ri = *loc.bins.entry((ci, b)).or_insert_with(|| { loc.rec.push((ci, b, vec![0; w])); loc.rec.len() - 1 });
            last = (rname.to_vec(), b, ri); ri
        };
        loc.reads += 1;
        for x in h.iter_mut() { *x = 0; }
        let (mut f, mut r, mut n, mut n_k) = (0u64, 0u64, 0usize, 0u64); kb.clear(); let cm = loc.cache.len() - 1;
        if seq != b"*" {
            for &c in seq {
                let c = LUT[c as usize];
                if c > 3 { n = 0; f = 0; r = 0; continue; }
                f = ((f << 2) | c as u64) & MASK; r = (r >> 2) | ((3 - c as u64) << 60); n += 1;
                if n >= K { n_k += 1; let x = f.min(r); let hx = hk(x); let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
                    if e.key == x { if e.sc != 0 { e.pend += 1; loc.hits += 1; let mut m = e.sc >> 25; while m != 0 { let j = m.trailing_zeros() as usize; h[j] += 1; m &= m - 1; } } }
                    else { prefetch(bloom.at(hx)); kb.push((x, hx, false)); } }
            }
        }
        loc.kmers += n_k;
        loc.cmiss += kb.len() as u64; for k in kb.iter_mut() { if bloom.has(k.1) { k.2 = true; prefetch(tab.at(k.1)); } }   // Bloom positives: prefetch their table slot
        for &(x, hx, bp) in kb.iter() {   // cache misses: evict, look up the table (Bloom positives only), fill
            let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
            if e.key != x { evict(e, tab); *e = match if bp { loc.bpos += 1; tab.find(x, hx) } else { None } { Some(i) => Ent { key: x, sc: i as u32 | (tab.s[i].meta >> 25) << 25, pend: 0 }, None => Ent { key: x, sc: 0, pend: 0 } }; }
            if e.sc != 0 { e.pend += 1; loc.hits += 1; let mut m = e.sc >> 25; while m != 0 { let j = m.trailing_zeros() as usize; h[j] += 1; m &= m - 1; } }
        }
        let rec = &mut loc.rec[ri].2;
        rec[0] += 1; if mq { rec[1] += 1; } if (fl >> 10) & 1 == 1 { rec[2] += 1; }
        let mut any = false;
        for j in 0..ncl { if h[j] > 0 { any = true; rec[3 + 3 * j] += h[j]; if mq { rec[4 + 3 * j] += h[j]; } rec[5 + 3 * j] += 1; } }
        if any { fa.push(b'>'); fa.extend_from_slice(rname); fa.push(b'|'); fa.extend_from_slice(b.to_string().as_bytes()); fa.push(b'|'); fa.extend_from_slice(mapq); fa.push(b'|'); fa.extend_from_slice(flag); fa.push(b'\n'); fa.extend_from_slice(seq); fa.push(b'\n'); }
    }
}

// Python repr of a float: shortest digits, fixed for exponents -4..15, else d.ddde-XX
fn pyf(x: f64) -> String {
    if x == 0.0 { return if x.is_sign_negative() { "-0.0".into() } else { "0.0".into() }; }
    if x.is_nan() { return "nan".into(); } if x.is_infinite() { return if x > 0.0 { "inf".into() } else { "-inf".into() }; }
    let e = format!("{x:e}"); let (m, ex) = e.split_once('e').unwrap(); let ex: i32 = ex.parse().unwrap();
    let neg = m.starts_with('-'); let d: String = m.trim_start_matches('-').chars().filter(|c| *c != '.').collect(); let n = d.len() as i32;
    let s = if (-4..16).contains(&ex) {
        if ex >= 0 { if n <= ex + 1 { format!("{d}{}.0", "0".repeat((ex + 1 - n) as usize)) } else { format!("{}.{}", &d[..(ex + 1) as usize], &d[(ex + 1) as usize..]) } }
        else { format!("0.{}{d}", "0".repeat((-ex - 1) as usize)) }
    } else { let mt = if n == 1 { d.clone() } else { format!("{}.{}", &d[..1], &d[1..]) }; format!("{mt}e{}{:02}", if ex < 0 { '-' } else { '+' }, ex.abs()) };
    if neg { format!("-{s}") } else { s }
}

// GC curve from ysinglecopy, applied to array classes
fn tm(a: &[u64], m: f64, lo: f64, hi: f64) -> f64 { let (mut s, mut n) = (0u64, 0u64); for &x in a { let y = x as f64; if lo * m <= y && y <= hi * m { s += x; n += 1; } } if n == 0 { f64::NAN } else { s as f64 / n as f64 } }
fn median(v: &[u64]) -> (f64, String) {   // statistics.median and its Python repr
    let mut s = v.to_vec(); let n = s.len(); let hi = *s.select_nth_unstable(n / 2).1;   // O(n) selection, same values as a sort
    if n % 2 == 1 { (hi as f64, hi.to_string()) } else { let lo = *s[..n / 2].iter().max().unwrap(); let m = (lo + hi) as f64 / 2.0; (m, pyf(m)) }
}
fn gc_sums(names: &[String], g: &[Vec<u8>], c: &[Vec<u64>]) -> Vec<f64> {
    let y = names.iter().position(|x| x == "ysinglecopy").unwrap_or_else(|| die("no ysinglecopy class"));
    let raw: Vec<f64> = c.iter().map(|v| v.iter().sum::<u64>() as f64).collect();
    let (m, _) = median(&c[y]); let al = tm(&c[y], m, 0.5, 1.5);
    if m < 10.0 || !(al > 0.0) { eprintln!("WARNING gc: ysinglecopy median count {m} < 10, no GC curve, sum_gc = sum"); return raw; }
    let mut by: Vec<Vec<u64>> = vec![vec![]; 32]; for (&q, &x) in g[y].iter().zip(&c[y]) { by[q as usize].push(x); }
    let d: Vec<(usize, f64)> = (0..32).filter(|&b| by[b].len() >= 1000).map(|b| (b, tm(&by[b], m, 0.25, 2.5) / al)).collect();
    if d.is_empty() || d.iter().any(|&(_, v)| !(v > 0.0)) { die("gc: bad curve"); }
    let r: Vec<f64> = (0..32).map(|b| { let mut best = d[0]; for &e in &d[1..] { if (e.0 as i64 - b as i64).abs() < (best.0 as i64 - b as i64).abs() { best = e; } } best.1 }).collect();
    let mut f = 0.0; for b in 0..32 { f += by[b].len() as f64 * r[b]; } f /= g[y].len() as f64;
    eprintln!("gc r(GC 10..20) {}", r[10..21].iter().map(|v| format!("{v:.3}")).collect::<Vec<_>>().join(" "));
    names.iter().enumerate().map(|(k, n)| if n != "ysinglecopy" && n != "chr2depth" { let mut s = 0.0; for (&q, &x) in g[k].iter().zip(&c[k]) { s += x as f64 / r[q as usize]; } s * f } else { raw[k] }).collect()
}

fn count(a: &[String]) {
    if a.len() != 6 { die("usage: samtools view ... | yruler count <markers.fa> <method.tsv or -> <name> <bin> <threads> <outdir>"); }
    let (fa_path, method, name, bin, thr, out) = (&a[0], &a[1], &a[2], a[3].parse::<u64>().unwrap_or_else(|_| die("bin")), a[4].parse::<usize>().unwrap_or_else(|_| die("threads")).max(1), &a[5]);
    let t0 = Instant::now();
    // reader: stdin read during the marker load (up to 64 chunks of 8 MB queued), chunks cut at line ends
    let (tx, rx) = mpsc::sync_channel::<(u64, Vec<u8>)>((thr * 2).max(64)); let rx = Arc::new(Mutex::new(rx));
    let (ptx, prx) = mpsc::channel::<Vec<u8>>();   // used chunk buffers back to the reader (no new pages per chunk)
    let reader = std::thread::spawn(move || {
        let mut carry: Vec<u8> = Vec::new(); let (mut ci, mut bytes) = (0u64, 0u64);
        loop {
            let mut buf = prx.try_recv().unwrap_or_else(|_| Vec::with_capacity(CHUNK + (1 << 17))); buf.clear(); buf.extend_from_slice(&carry); carry.clear();
            let start = buf.len(); buf.reserve(CHUNK); let bp = unsafe { buf.as_mut_ptr().add(start) };   // read(2) straight into spare capacity: no zero fill
            let mut got = 0; while got < CHUNK { let k = unsafe { libc::read(0, bp.add(got) as *mut libc::c_void, CHUNK - got) }; if k == 0 { break; } if k < 0 { let e = io::Error::last_os_error(); if e.kind() == io::ErrorKind::Interrupted { continue; } die(&format!("stdin: {e}")); } got += k as usize; }
            unsafe { buf.set_len(start + got); } bytes += got as u64;
            if got == 0 { if !buf.is_empty() { tx.send((ci, buf)).unwrap(); } break; }
            if let Some(p) = buf.iter().rposition(|&b| b == b'\n') { carry.extend_from_slice(&buf[p + 1..]); buf.truncate(p + 1); }
            tx.send((ci, buf)).unwrap(); ci += 1;
        }
        bytes
    });
    let mut mk = load(fa_path, thr); let ncl = mk.names.len(); let n = mk.key.len();
    let mt: HashMap<String, String> = if method == "-" { String::new() } else { std::fs::read_to_string(method).unwrap_or_else(|e| die(&format!("{method}: {e}"))) }.lines().filter(|l| !l.starts_with('#')).filter_map(|l| { let mut s = l.split('\t'); Some((s.next()?.to_string(), s.next()?.to_string())) }).collect();
    let per: Vec<usize> = (0..ncl).map(|j| mk.ocls.iter().filter(|&&c| c == j).count()).collect();
    for (j, c) in mk.names.iter().enumerate() { if let Some(e) = mt.get(&format!("n_{c}")) { if e.parse::<usize>().ok() != Some(per[j]) { die(&format!("{c} k-mers {} vs method {e}", per[j])); } } }
    let tab = std::mem::replace(&mut mk.tab, Table { s: Vec::new(), mask: 0 });
    eprintln!("markers lines {} distinct {} classes {:?} per class {:?}; table {} MB; load {:.1} s", mk.order.len(), n, mk.names, per, (tab.s.len() * 16) >> 20, t0.elapsed().as_secs_f64());
    let bpk: usize = std::env::var("YRULER_BLOOM_BITS").ok().and_then(|v| v.parse().ok()).unwrap_or(16);   // Bloom bits per key
    let bloom = { let mut b = Bloom::new(n, bpk.max(1)); for &x in &mk.key { b.insert(hk(x)); } Arc::new(b) };
    eprintln!("bloom {} MB ({} bits/key)", bloom.b.len() * 8 >> 20, bpk.max(1));
    if tab.s.len() > (IDX as usize) + 1 { die("table above 2^25 slots: cache entry slot field is 25 bits"); }
    let tab = Arc::new(tab);
    let cb: u32 = std::env::var("YRULER_CACHE_BITS").ok().and_then(|v| v.parse().ok()).unwrap_or(15);   // per-thread cache: 2^cb entries
    let (ftx, frx) = mpsc::sync_channel::<(u64, Vec<u8>)>(thr * 4);
    std::fs::create_dir_all(out).unwrap_or_else(|e| die(&format!("{out}: {e}")));
    let any_path = format!("{out}/any.fa");
    let writer = std::thread::spawn(move || { let mut w = BufWriter::with_capacity(1 << 22, std::fs::File::create(&any_path).unwrap_or_else(|e| die(&format!("{any_path}: {e}")))); let mut next = 0u64; let mut held: std::collections::BTreeMap<u64, Vec<u8>> = std::collections::BTreeMap::new();   // marker reads in input order
        let unordered = std::env::var("YRULER_UNORDERED").is_ok();   // write blocks as they finish
        for (i, b) in frx { if unordered { w.write_all(&b).unwrap_or_else(|e| die(&format!("write any.fa: {e}"))); continue; } held.insert(i, b); while let Some(b) = held.remove(&next) { w.write_all(&b).unwrap_or_else(|e| die(&format!("write any.fa: {e}"))); next += 1; } }
        if !held.is_empty() { die("marker reads: chunks missing"); } w.flush().unwrap(); });
    let workers: Vec<_> = (0..thr).map(|_| {
        let (rx, ftx, ptx, tab, bloom) = (rx.clone(), ftx.clone(), ptx.clone(), tab.clone(), bloom.clone());
        std::thread::spawn(move || {
            let mut loc = Local { contigs: HashMap::new(), cname: vec![], bins: HashMap::new(), rec: vec![], reads: 0, kmers: 0, hits: 0, cmiss: 0, bpos: 0, cache: vec![Ent { key: EMPTY, sc: 0, pend: 0 }; 1 << cb] }; let mut kb: Vec<(u64, u64, bool)> = Vec::with_capacity(1024);
            loop {
                let c = { let r = rx.lock().unwrap(); r.recv() };
                let Ok((ci, c)) = c else { break };
                let mut fa = Vec::with_capacity(1 << 16);
                work(&c, &mut loc, &bloom, &tab, ncl, bin, &mut fa, &mut kb);
                ftx.send((ci, fa)).unwrap(); let _ = ptx.send(c);   // fails only after the reader ended: buffer then freed
            }
            for e in &loc.cache { evict(e, &tab); }
            loc
        })
    }).collect();
    drop(ftx);
    let bytes_in = reader.join().unwrap();
    let locs: Vec<Local> = workers.into_iter().map(|w| w.join().unwrap()).collect();
    writer.join().unwrap();
    let t_count = t0.elapsed().as_secs_f64();
    let tp = Instant::now(); let lap = |w: &str| eprintln!("post {w} at {:.2} s", tp.elapsed().as_secs_f64());
    // merge bins
    let w = 3 + 3 * ncl; let mut bins: HashMap<(Vec<u8>, u64), Vec<u64>> = HashMap::new();
    let (mut reads, mut kmers, mut hits, mut cmiss, mut bpos) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for l in &locs {
        reads += l.reads; kmers += l.kmers; cmiss += l.cmiss; bpos += l.bpos; hits += l.hits;
        for (ci, b, v) in &l.rec { let e = bins.entry((l.cname[*ci as usize].clone(), *b)).or_insert_with(|| vec![0; w]); for (x, y) in e.iter_mut().zip(v) { *x += y; } }
    }
    if reads == 0 { die("0 reads"); }
    eprintln!("reads {reads} bins {} bytes {} k-mers {kmers} cache misses {cmiss} ({:.2}%) bloom positives {bpos} hits {hits}; counting {:.1} s", bins.len(), bytes_in, 100.0 * cmiss as f64 / kmers.max(1) as f64, t_count);
    lap("bins merged");
    let mut cnt: Vec<u64> = vec![0; n]; for sl in tab.s.iter() { if sl.key != EMPTY { cnt[(sl.meta & IDX) as usize] = sl.cnt.load(Relaxed) as u64; } }
    // per class lists in marker order
    let mut cv: Vec<Vec<u64>> = vec![vec![]; ncl]; let mut gv: Vec<Vec<u8>> = vec![vec![]; ncl];
    for (k, &i) in mk.order.iter().enumerate() { let c = mk.ocls[k]; cv[c].push(cnt[i as usize]); gv[c].push(mk.gc[k]); }
    // check: bins hits = counts per class
    for j in 0..ncl { let s: u64 = cv[j].iter().sum(); let bh: u64 = bins.values().map(|v| v[3 + 3 * j]).sum(); eprintln!("check {} counts {s} bins {bh}", mk.names[j]); if s != bh { die(&format!("{} bins hits {bh} vs counts {s}", mk.names[j])); } }
    let wr = |p: &str| BufWriter::with_capacity(1 << 22, std::fs::File::create(format!("{out}/{p}")).unwrap_or_else(|e| die(&format!("{p}: {e}"))));
    let t1 = Instant::now();
    { // counts.tsv: text built in parallel blocks, written in order
        let per = n.div_ceil(thr).max(1); let (key, cnt) = (&mk.key, &cnt);
        let parts: Vec<Vec<u8>> = std::thread::scope(|sc| { let hs: Vec<_> = (0..n.div_ceil(per)).map(|t| sc.spawn(move || { let mut b = Vec::with_capacity(per * 40); for i in t * per..((t + 1) * per).min(n) { b.extend_from_slice(&decode(key[i])); b.push(b'\t'); b.extend_from_slice(cnt[i].to_string().as_bytes()); b.push(b'\n'); } b })).collect(); hs.into_iter().map(|h| h.join().unwrap()).collect() });
        let mut o = wr("counts.tsv"); for b in &parts { o.write_all(b).unwrap_or_else(|e| die(&format!("write counts.tsv: {e}"))); } o.flush().unwrap();
    }
    eprintln!("write counts {:.2} s", t1.elapsed().as_secs_f64());
    {
        let mut o = wr("bins.tsv"); let mut keys: Vec<_> = bins.keys().cloned().collect(); keys.sort();
        write!(o, "contig\tstart\treads\treads_mq20\tduplicates").unwrap();
        for c in &mk.names { write!(o, "\t{c}_hits\t{c}_hits_mq20\t{c}_reads").unwrap(); } writeln!(o).unwrap();
        for k in keys { let v = &bins[&k]; o.write_all(&k.0).unwrap(); write!(o, "\t{}", k.1 * bin).unwrap(); for x in v { write!(o, "\t{x}").unwrap(); } writeln!(o).unwrap(); }
        o.flush().unwrap();
    }
    let need = ["ysinglecopy", "chr2depth"]; for c in need { if !mk.names.iter().any(|x| x == c) { die(&format!("missing class {c}")); } }
    let ix = |c: &str| mk.names.iter().position(|x| x == c).unwrap();
    let xd = cv[ix("ysinglecopy")].iter().sum::<u64>() as f64 / cv[ix("ysinglecopy")].len() as f64;
    let ad = cv[ix("chr2depth")].iter().sum::<u64>() as f64 / cv[ix("chr2depth")].len() as f64;
    if ad == 0.0 { die("chr2depth mean 0"); }
    lap("counts and bins written");
    let sg = gc_sums(&mk.names, &gv, &cv); lap("gc");
    {
        let mut o = wr("summary.tsv");
        writeln!(o, "#y_dose={}\tysinglecopy_depth={}\tchr2depth_mean={}", pyf(2.0 * xd / ad), pyf(xd), pyf(ad)).unwrap();
        writeln!(o, "sample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc").unwrap();
        for (j, c) in mk.names.iter().enumerate() { let v = &cv[j]; let (_, ms) = median(v); let z = v.iter().filter(|&&x| x == 0).count(); writeln!(o, "{name}\t{reads}\t{c}\t{}\t{}\t{ms}\t{:.4}\t{:.1}", v.len(), v.iter().sum::<u64>(), z as f64 / v.len() as f64, sg[j]).unwrap(); }
        o.flush().unwrap();
    }
    eprintln!("y_dose {:?} ysinglecopy_depth {:?} chr2depth_mean {:?}; total {:.1} s", 2.0 * xd / ad, xd, ad, t0.elapsed().as_secs_f64());
}

// estimate: length = K * sum_gc(class) / sum(ysinglecopy)
fn estimate(a: &[String]) {
    if a.len() != 2 { die("usage: yruler estimate <summary.tsv> <model.tsv>"); }
    let rd = |p: &str| std::fs::read_to_string(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let (sm, md) = (rd(&a[0]), rd(&a[1]));
    let mut dose = String::new(); let mut sg: Vec<(String, String, f64)> = vec![]; let mut ysum = f64::NAN;
    for l in sm.lines() {
        if let Some(h) = l.strip_prefix('#') { dose = h.split('\t').find_map(|f| f.strip_prefix("y_dose=")).unwrap_or("").to_string(); continue; }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 8 || f[0] == "sample" { continue; }
        let v: f64 = f[7].parse().unwrap_or_else(|_| die(&format!("sum_gc {}", f[7])));
        if f[2] == "ysinglecopy" { ysum = f[4].parse().unwrap_or_else(|_| die("ysinglecopy sum")); }
        sg.push((f[0].to_string(), f[2].to_string(), v));
    }
    if !(ysum > 0.0) { die("no ysinglecopy sum > 0 in summary"); }
    let mut k: HashMap<String, f64> = HashMap::new();
    for l in md.lines().filter(|l| !l.starts_with('#')) { let f: Vec<&str> = l.split('\t').collect(); if f.len() >= 2 && f[0] != "class" { k.insert(f[0].to_string(), f[1].parse().unwrap_or_else(|_| die(&format!("model K {}", f[1])))); } }
    if k.is_empty() { die("model has no class rows"); }
    println!("sample\tclass\tlength_bp\tmarker_sum_gc\tysinglecopy_sum\ty_dose");
    let mut n = 0;
    for (s, c, v) in &sg { if let Some(kk) = k.get(c) { println!("{s}\t{c}\t{:.0}\t{v:.1}\t{ysum:.0}\t{dose}", kk * v / ysum); n += 1; } }
    for c in k.keys() { if !sg.iter().any(|t| &t.1 == c) { die(&format!("model class {c} missing in summary")); } }
    eprintln!("estimate classes {n} y_dose {dose}");
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("count") => count(&a[2..]),
        Some("estimate") => estimate(&a[2..]),
        Some("--version") => println!("yruler {}", env!("CARGO_PKG_VERSION")),
        _ => die("usage: yruler count|estimate ... (see README)"),
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn masked_superstring() {   // uppercase starts a k-mer
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGTa";
        let mut ms = b">c1\n".to_vec(); ms.extend(s.iter().enumerate().map(|(i, b)| if i < 2 { b.to_ascii_uppercase() } else { b.to_ascii_lowercase() })); ms.push(b'\n');
        let e = super::expand_ms(&ms); let l: Vec<&[u8]> = e.split(|&b| b == b'\n').filter(|l| !l.is_empty()).collect();
        assert_eq!(l[0], b">c1"); assert_eq!(l[1], &s[0..31]); assert_eq!(l[2], b"CGTACGTACGTACGTACGTACGTACGTACGT"); assert_eq!(l.len(), 3);
    }
    #[test]
    fn median_select() {   // O(n) median = sort median
        let mut x = 7u64; for n in 1..400 { let v: Vec<u64> = (0..n).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (x >> 59) % 13 }).collect();
            let mut s = v.clone(); s.sort(); let m = if n % 2 == 1 { s[n / 2] as f64 } else { (s[n / 2 - 1] + s[n / 2]) as f64 / 2.0 }; assert_eq!(super::median(&v).0, m); }
    }
    #[test]
    fn python_repr() {   // values from python3 repr()
        for (x, s) in [(4.517420114170894e-05, "4.517420114170894e-05"), (0.0001756600924138287, "0.0001756600924138287"), (1.0095530193019127, "1.0095530193019127"), (27.94310690421051, "27.94310690421051"), (14.0, "14.0"), (1e16, "1e+16"), (123456789012345678.0, "1.2345678901234568e+17"), (0.0001, "0.0001"), (2.5e-10, "2.5e-10"), (14.5, "14.5")] {
            assert_eq!(super::pyf(x), s);
        }
    }
}

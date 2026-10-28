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
const CHUNK: usize = 2 << 20;   // read chunk bytes (2 MB: same speed as 8 to 64, least RSS); a longer line grows its chunk
fn chunk() -> usize { CHUNK }

static LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'C' as usize] = 1; t[b'G' as usize] = 2; t[b'T' as usize] = 3;
    t[b'a' as usize] = 0; t[b'c' as usize] = 1; t[b'g' as usize] = 2; t[b't' as usize] = 3;
    t
};


#[inline(always)]
fn hk(x: u64) -> u64 { let m = x.wrapping_mul(0x9e3779b97f4a7c15); m ^ (m >> 29) }   // k-mer hash: one multiply, high bits folded down

fn die(m: &str) -> ! { eprintln!("ERROR yruler: {m}"); std::process::exit(1) }
// resident and peak memory from /proc (Linux), for the stage log
fn rss() -> String { std::fs::read_to_string("/proc/self/status").map_or("rss n/a".into(), |t| t.lines().filter(|l| l.starts_with("VmRSS") || l.starts_with("VmHWM")).map(|l| l.split_whitespace().take(2).collect::<Vec<_>>().join(" ")).collect::<Vec<_>>().join(" ")) }

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

// open-addressing table: 16-byte slot (key, count, marker index), class bits per slot (up to 16 classes)
#[repr(C, align(16))]
struct Slot { key: u64, cnt: AtomicU32, idx: u32 }
fn hugepage<T>(v: &Vec<T>) {   // ask the kernel for transparent huge pages before the first touch (fewer TLB misses)
    let (p, n) = (v.as_ptr() as usize, v.capacity() * std::mem::size_of::<T>()); let a = (p + (1 << 21) - 1) & !((1 << 21) - 1);
    if n > (4 << 21) { unsafe { libc::madvise(a as *mut libc::c_void, (p + n - a) & !((1 << 21) - 1), libc::MADV_HUGEPAGE); } }
}
struct Table { s: Vec<Slot>, cls: Vec<u16>, mask: usize }
impl Table {
    fn new(n: usize, thr: usize) -> Self {   // slots written in parallel (first touch spread over threads)
        let s = (n * 2).next_power_of_two(); let mut v: Vec<Slot> = Vec::with_capacity(s); hugepage(&v);
        let (p, per) = (v.as_mut_ptr() as usize, s.div_ceil(thr));
        std::thread::scope(|sc| for t in 0..thr { sc.spawn(move || { let p = p as *mut Slot; for i in t * per..((t + 1) * per).min(s) { unsafe { p.add(i).write(Slot { key: EMPTY, cnt: AtomicU32::new(0), idx: 0 }); } } }); });
        unsafe { v.set_len(s); } Table { s: v, cls: vec![0; s], mask: s - 1 }
    }
    // returns the marker index of key, inserting it with index val; class bit c is ORed in
    fn insert(&mut self, key: u64, h: u64, val: u32, c: usize) -> (u32, bool) {
        let mut i = (h >> 7) as usize & self.mask;
        loop {
            let sl = &mut self.s[i];
            if sl.key == EMPTY { sl.key = key; sl.idx = val; self.cls[i] |= 1 << c; return (val, true); }
            if sl.key == key { self.cls[i] |= 1 << c; return (sl.idx, false); }
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

struct Markers { key: Vec<u64>, order: Vec<u32>, ocls: Vec<u8>, gc: Vec<u8>, names: Vec<String>, tab: Table }
// class name of a header line: text before the first '_' (trailing whitespace trimmed first)
fn hname(h: &[u8]) -> &[u8] { let mut e = h.len(); while e > 0 && h[e - 1].is_ascii_whitespace() { e -= 1; } let h = &h[..e]; &h[..h.iter().position(|&b| b == b'_').unwrap_or(h.len())] }
// masked superstring records (class name, sequence): uppercase base = a k-mer starts there
fn ms_records(raw: &[u8]) -> Vec<(&[u8], Vec<u8>)> {
    let mut r: Vec<(&[u8], Vec<u8>)> = vec![];
    for l in raw.split(|&b| b == b'\n') { let l = l.strip_suffix(b"\r").unwrap_or(l); if l.first() == Some(&b'>') { r.push((hname(&l[1..]), vec![])); } else if let Some(x) = r.last_mut() { x.1.extend_from_slice(l); } else if !l.is_empty() { die("sequence before the first header"); } }
    r
}
type Part<'a> = (Vec<(u64, u8, u8)>, Vec<&'a [u8]>);   // (key, G+C, local header index or 255), header names
fn gcn(w: &[u8]) -> u8 { w.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count() as u8 }
// k-mers in file order, parsed in parallel pieces; first occurrence = marker index
fn load(path: &str, thr: usize) -> Markers {
    let raw = std::fs::read(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    let raw = if path.ends_with(".zst") { zstd::stream::decode_all(&raw[..]).unwrap_or_else(|e| die(&format!("{path}: zstd {e}"))) } else { raw };
    // one k-mer per line, else masked superstring
    let ms = raw.split(|&b| b == b'\n').find(|l| !l.is_empty() && l[0] != b'>').map_or(false, |l| { let l = l.strip_suffix(b"\r").unwrap_or(l); l.len() != K || l.iter().any(|b| b.is_ascii_lowercase()) });
    let recs = if ms { ms_records(&raw) } else { vec![] };
    let parts: Vec<Part> = if ms {   // record pieces of about starts/thr k-mer starts
        let ns: usize = recs.iter().map(|r| (r.1.len() + 1).saturating_sub(K)).sum(); let per = (ns / thr).max(1 << 16);
        let tasks: Vec<(usize, usize, usize)> = recs.iter().enumerate().flat_map(|(ri, r)| { let n = (r.1.len() + 1).saturating_sub(K); (0..n).step_by(per).map(move |a| (ri, a, (a + per).min(n))) }).collect();
        let p: Vec<Part> = std::thread::scope(|sc| { let recs = &recs; let hs: Vec<_> = tasks.iter().map(|&(ri, a, b)| sc.spawn(move || { let s = &recs[ri].1; let mut v = Vec::with_capacity((b - a) / 2);
            for i in a..b { if s[i].is_ascii_uppercase() { let w = &s[i..i + K]; v.push((canon(w).unwrap_or_else(|| die("masked superstring k-mer with non-ACGT")), gcn(w), 0)); } }
            (v, vec![recs[ri].0]) })).collect(); hs.into_iter().map(|h| h.join().unwrap()).collect() });
        eprintln!("markers: masked superstring {} bytes, {} records, {} k-mer starts, {} pieces", raw.len(), recs.len(), ns, tasks.len()); p
    } else {   // file cut at line ends into thr pieces
        let n = raw.len();
        let mut cuts = vec![0usize]; for t in 1..thr { let mut p = (n * t / thr).max(*cuts.last().unwrap()); while p < n && p > 0 && raw[p - 1] != b'\n' { p += 1; } cuts.push(p); } cuts.push(n);
        std::thread::scope(|sc| {
            let hs: Vec<_> = cuts.windows(2).map(|w| { let raw = &raw[w[0]..w[1]]; sc.spawn(move || {
                let (mut v, mut nm, mut cur): (Vec<(u64, u8, u8)>, Vec<&[u8]>, u8) = (Vec::with_capacity(raw.len() / 64 + 16), vec![], 255);
                for l in raw.split(|&b| b == b'\n') {
                    if l.is_empty() { continue; }
                    if l[0] == b'>' { let c = hname(&l[1..]); if nm.last() != Some(&c) { if nm.len() == 254 { die("more than 254 header names in one piece"); } nm.push(c); } cur = (nm.len() - 1) as u8; continue; }
                    let l = l.strip_suffix(b"\r").unwrap_or(l);
                    v.push((canon(l).unwrap_or_else(|| die(&format!("marker line of length {} or non-ACGT", l.len()))), gcn(l), cur));
                }
                (v, nm) }) }).collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        })
    };
    let mut names: Vec<String> = vec![]; let mut cur = usize::MAX; let tot: usize = parts.iter().map(|p| p.0.len()).sum();
    if tot >= u32::MAX as usize { die("more than 2^32 - 1 markers"); }
    let mut tab = Table::new(tot, thr); let (mut key, mut order, mut ocls, mut gc) = (vec![], Vec::with_capacity(tot), Vec::with_capacity(tot), Vec::with_capacity(tot));
    for (v, nm) in &parts {
        let gl: Vec<usize> = nm.iter().map(|c| { let c = String::from_utf8_lossy(c).to_string(); names.iter().position(|x| *x == c).unwrap_or_else(|| { names.push(c); names.len() - 1 }) }).collect();
        if names.len() > 16 { die("more than 16 classes"); }
        for (k, &(x, g, lc)) in v.iter().enumerate() {
            if lc != 255 { cur = gl[lc as usize]; } else if cur == usize::MAX { die("sequence before the first header"); }
            if let Some(nx) = v.get(k + 16) { prefetch(tab.at(hk(nx.0))); }
            let (i, new) = tab.insert(x, hk(x), key.len() as u32, cur); if new { key.push(x); }
            order.push(i); ocls.push(cur as u8); gc.push(g);
        }
        if let Some(&c) = gl.last() { cur = c; }
    }
    Markers { key, order, ocls, gc, names, tab }
}

// per-thread bin record (rk key, rv flat values): reads, reads MAPQ >= 20, duplicates, then per class hits, hits MAPQ >= 20, reads with a hit
struct Local { contigs: HashMap<Vec<u8>, u32>, cname: Vec<Vec<u8>>, bins: HashMap<(u32, u64), usize>, rk: Vec<(u32, u64)>, rv: Vec<u64>, reads: u64, kmers: u64, hits: u64, cmiss: u64, bpos: u64, cache: Vec<Ent> }
// per-thread cache of recent k-mers (sorted input: a k-mer comes back in the next reads); counts added to the table on eviction
#[derive(Clone, Copy)]
struct Ent { key: u64, sc: u32, pend: u16, cls: u16 }   // 16 bytes; sc = table slot + 1, 0 = not a marker
fn evict(e: &Ent, tab: &Table) { if e.sc != 0 && e.pend > 0 { tab.s[e.sc as usize - 1].cnt.fetch_add(e.pend as u32, Relaxed); } }
#[inline(always)]
fn hit(e: &mut Ent, tab: &Table, h: &mut [u64]) {   // count, flush before u16 overflow
    if e.pend == u16::MAX { tab.s[e.sc as usize - 1].cnt.fetch_add(e.pend as u32, Relaxed); e.pend = 0; }
    e.pend += 1; let mut m = e.cls; while m != 0 { h[m.trailing_zeros() as usize] += 1; m &= m - 1; }
}

#[inline(always)]
fn mc(s: &[u8], c: u8) -> Option<usize> { let p = unsafe { libc::memchr(s.as_ptr() as *const libc::c_void, c as i32, s.len()) }; if p.is_null() { None } else { Some(p as usize - s.as_ptr() as usize) } }   // glibc SIMD memchr
#[inline(always)]
fn field<'a>(l: &'a [u8], q: &mut usize) -> &'a [u8] {   // next tab-separated field from q (last field: to line end)
    if *q > l.len() { die("SAM line with too few fields"); }
    let r = &l[*q..]; let e = mc(r, b'\t').unwrap_or(r.len()); *q += e + 1; &r[..e]
}
fn num(s: &[u8]) -> u64 { let mut v = 0u64; for &b in s { if !b.is_ascii_digit() { die("non-numeric SAM field"); } v = v * 10 + (b - b'0') as u64; } v }

fn itoa(b: &mut Vec<u8>, mut v: u64) { let mut t = [0u8; 20]; let mut i = 20; loop { i -= 1; t[i] = b'0' + (v % 10) as u8; v /= 10; if v == 0 { break; } } b.extend_from_slice(&t[i..]); }
struct Ctx<'a> { bloom: &'a Bloom, tab: &'a Table, ncl: usize, bin: u64, h: Vec<u64>, last: (Vec<u8>, u64, usize) }
// one read: k-mer hits, bin record, marker read
fn read1(c: &mut Ctx, loc: &mut Local, fa: &mut Vec<u8>, kb: &mut Vec<(u64, u64, bool)>, rname: &[u8], pos: u64, mapq: u64, fl: u64, seq: &[u8]) {
    let (w, b, mq) = (3 + 3 * c.ncl, pos / c.bin, mapq >= 20);
    let ri = if c.last.1 == b && c.last.0 == rname { c.last.2 } else {
        let ci = match loc.contigs.get(rname) { Some(&x) => x, None => { let x = loc.cname.len() as u32; loc.contigs.insert(rname.to_vec(), x); loc.cname.push(rname.to_vec()); x } };
        let ri = *loc.bins.entry((ci, b)).or_insert_with(|| { loc.rk.push((ci, b)); loc.rv.resize(loc.rv.len() + w, 0); loc.rk.len() - 1 });
        c.last = (rname.to_vec(), b, ri); ri
    };
    loc.reads += 1; c.h.iter_mut().for_each(|x| *x = 0); let h = &mut c.h;
    let (mut f, mut r, mut n, mut n_k) = (0u64, 0u64, 0usize, 0u64); kb.clear(); let cm = loc.cache.len() - 1;
    if seq != b"*" {
        for &ch in seq {
            let ch = LUT[ch as usize];
            if ch > 3 { n = 0; f = 0; r = 0; continue; }
            f = ((f << 2) | ch as u64) & MASK; r = (r >> 2) | ((3 - ch as u64) << 60); n += 1;
            if n >= K { n_k += 1; let x = f.min(r); let hx = hk(x); let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
                if e.key == x { if e.sc != 0 { loc.hits += 1; hit(e, c.tab, h); } }
                else { prefetch(c.bloom.at(hx)); kb.push((x, hx, false)); } }
        }
    }
    loc.kmers += n_k;
    loc.cmiss += kb.len() as u64; for k in kb.iter_mut() { if c.bloom.has(k.1) { k.2 = true; prefetch(c.tab.at(k.1)); } }   // Bloom positives: prefetch slot
    for &(x, hx, bp) in kb.iter() {   // cache misses: evict, look up, fill
        let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
        if e.key != x { evict(e, c.tab); *e = match if bp { loc.bpos += 1; c.tab.find(x, hx) } else { None } { Some(i) => Ent { key: x, sc: i as u32 + 1, pend: 0, cls: c.tab.cls[i] }, None => Ent { key: x, sc: 0, pend: 0, cls: 0 } }; }
        if e.sc != 0 { loc.hits += 1; hit(e, c.tab, h); }
    }
    let rec = &mut loc.rv[ri * w..(ri + 1) * w];
    rec[0] += 1; if mq { rec[1] += 1; } if (fl >> 10) & 1 == 1 { rec[2] += 1; }
    let mut any = false;
    for j in 0..c.ncl { if h[j] > 0 { any = true; rec[3 + 3 * j] += h[j]; if mq { rec[4 + 3 * j] += h[j]; } rec[5 + 3 * j] += 1; } }
    if any { fa.push(b'>'); fa.extend_from_slice(rname); fa.push(b'|'); itoa(fa, b); fa.push(b'|'); itoa(fa, mapq); fa.push(b'|'); itoa(fa, fl); fa.push(b'\n'); fa.extend_from_slice(seq); fa.push(b'\n'); }
}
// SAM text chunk
fn work(chunk: &[u8], loc: &mut Local, c: &mut Ctx, fa: &mut Vec<u8>, kb: &mut Vec<(u64, u64, bool)>) {
    let mut p0 = 0;
    while p0 < chunk.len() {
        let e0 = mc(&chunk[p0..], b'\n').map_or(chunk.len(), |i| p0 + i); let line = &chunk[p0..e0]; p0 = e0 + 1;
        if line.is_empty() || line[0] == b'@' { continue; }
        let mut q = 0; let it = &mut q;
        field(line, it); let fl = num(field(line, it)); let rname = field(line, it); let pos = num(field(line, it)); let mapq = num(field(line, it));
        for _ in 0..4 { field(line, it); } let seq = field(line, it);
        read1(c, loc, fa, kb, rname, pos, mapq, fl, seq);
    }
}
// binary chunk (first byte 0xFF): flag u16, tid i32, pos i64, mapq u8, l_seq u32, 4-bit seq
fn work_bin(ch: &[u8], loc: &mut Local, c: &mut Ctx, fa: &mut Vec<u8>, kb: &mut Vec<(u64, u64, bool)>, sq: &mut Vec<u8>) {
    const NT: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let names = NAMES.get().unwrap_or_else(|| die("contig names not set")); let mut p = 1;
    while p < ch.len() {
        let fl = u16::from_le_bytes([ch[p], ch[p + 1]]) as u64; let tid = i32::from_le_bytes(ch[p + 2..p + 6].try_into().unwrap()); let pos = (i64::from_le_bytes(ch[p + 6..p + 14].try_into().unwrap()) + 1).max(0) as u64;
        let mq = ch[p + 14] as u64; let n = u32::from_le_bytes(ch[p + 15..p + 19].try_into().unwrap()) as usize; p += 19; let e = &ch[p..p + (n + 1) / 2]; p += (n + 1) / 2;
        sq.clear(); if n == 0 { sq.push(b'*'); } else { sq.extend((0..n).map(|i| NT[((e[i >> 1] >> (4 * (1 - (i & 1)))) & 15) as usize])); }
        read1(c, loc, fa, kb, if tid < 0 { b"*" } else { &names[tid as usize] }, pos, mq, fl, sq);
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

static NAMES: std::sync::OnceLock<Vec<Vec<u8>>> = std::sync::OnceLock::new();
enum Src { Stdin, Hts { input: String, reference: Option<String>, bed: Option<String> } }
type Tx = mpsc::SyncSender<(u64, Vec<u8>)>;
// reads into SAM-like chunks cut at line ends
fn spawn_reader(src: Src, tx: Tx, prx: mpsc::Receiver<Vec<u8>>, thr: usize) -> std::thread::JoinHandle<u64> {
    std::thread::Builder::new().name("yr-read".into()).spawn(move || match src { Src::Stdin => read_stdin(tx, prx), Src::Hts { input, reference, bed } => read_hts(&input, reference.as_deref(), bed.as_deref(), thr, tx, prx) }).unwrap()
}
fn read_stdin(tx: Tx, prx: mpsc::Receiver<Vec<u8>>) -> u64 {
    let mut carry: Vec<u8> = Vec::new(); let (mut ci, mut bytes) = (0u64, 0u64);
    loop {
        let mut buf = prx.try_recv().unwrap_or_else(|_| Vec::with_capacity(chunk() + (1 << 17))); buf.clear(); buf.extend_from_slice(&carry); carry.clear();
        let start = buf.len(); buf.reserve(chunk()); let bp = unsafe { buf.as_mut_ptr().add(start) };   // read(2) into spare capacity
        let mut got = 0; while got < chunk() { let k = unsafe { libc::read(0, bp.add(got) as *mut libc::c_void, chunk() - got) }; if k == 0 { break; } if k < 0 { let e = io::Error::last_os_error(); if e.kind() == io::ErrorKind::Interrupted { continue; } die(&format!("stdin: {e}")); } got += k as usize; }
        unsafe { buf.set_len(start + got); } bytes += got as u64;
        if got == 0 { if !buf.is_empty() { tx.send((ci, buf)).unwrap(); } break; }
        match buf.iter().rposition(|&b| b == b'\n') { Some(p) => { carry.extend_from_slice(&buf[p + 1..]); buf.truncate(p + 1); } None => { carry = buf; continue; } }   // no line end: line longer than the chunk, read on
        tx.send((ci, buf)).unwrap(); ci += 1;
    }
    bytes
}
// BED regions merged per contig, in header order
fn regions(bed: &str, names: &[Vec<u8>]) -> Vec<(u32, i64, i64)> {
    let mut by: HashMap<Vec<u8>, Vec<(i64, i64)>> = HashMap::new();
    for l in std::fs::read_to_string(bed).unwrap_or_else(|e| die(&format!("{bed}: {e}"))).lines() {
        let f: Vec<&str> = l.split('\t').collect(); if f.len() < 3 || f[0].starts_with('#') || f[0] == "*" { continue; }
        let (a, b) = (f[1].parse::<i64>().unwrap_or_else(|_| die(&format!("bed start {}", f[1]))), f[2].parse::<i64>().unwrap_or_else(|_| die(&format!("bed end {}", f[2]))));
        by.entry(f[0].as_bytes().to_vec()).or_default().push((a, b));
    }
    let mut out = vec![];
    for (t, n) in names.iter().enumerate() { if let Some(v) = by.get_mut(n) { v.sort(); let mut m: Vec<(i64, i64)> = vec![]; for &(a, b) in v.iter() { match m.last_mut() { Some(l) if a <= l.1 => l.1 = l.1.max(b), _ => m.push((a, b)) } } for (a, b) in m { out.push((t as u32, a, b)); } } }
    let found: usize = by.keys().filter(|k| names.contains(k)).count(); if found < by.len() { eprintln!("WARNING bed: {} of {} contigs not in the input header", by.len() - found, by.len()); }
    out
}
// CRAM/BAM: shards (region, previous end on contig) read in parallel, each read once
fn read_hts(input: &str, reference: Option<&str>, bed: Option<&str>, thr: usize, tx: Tx, prx: mpsc::Receiver<Vec<u8>>) -> u64 {
    use rust_htslib::bam::{self, Read};
    use std::sync::atomic::{AtomicU64, AtomicUsize};
    let open = |t: usize| { let mut r = if input.contains("://") { bam::IndexedReader::from_url(&input.parse().unwrap_or_else(|e| die(&format!("{input}: {e}")))) } else { bam::IndexedReader::from_path(input) }.unwrap_or_else(|e| die(&format!("{input}: {e}")));
        if let Some(f) = reference { r.set_reference(f).unwrap_or_else(|e| die(&format!("reference {f}: {e}"))); }
        r.set_threads(t.max(1)).unwrap_or_else(|e| die(&format!("threads: {e}")));
        if let Err(e) = r.set_cram_options(rust_htslib::htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, 0x2 | 0x4 | 0x8 | 0x10 | 0x200) { eprintln!("WARNING cram required fields: {e}"); }
        r };
    let h = open(1); let names: Vec<Vec<u8>> = h.header().target_names().iter().map(|n| n.to_vec()).collect(); let _ = NAMES.set(names.clone());
    let rg: Vec<(u32, i64, i64)> = match bed { Some(b) => regions(b, &names),
        None => (0..names.len() as u32).flat_map(|t| { let l = h.header().target_len(t).unwrap_or(0) as i64; (0..l.max(1)).step_by(10_000_000).map(move |a| (t, a, (a + 10_000_000).min(l))) }).collect() };
    let mut sh: Vec<Option<(u32, i64, i64, i64)>> = rg.iter().enumerate().map(|(i, &(t, a, e))| Some((t, a, e, if i > 0 && rg[i - 1].0 == t { rg[i - 1].2 } else { i64::MIN }))).collect();
    sh.push(None);   // unmapped
    let nr = (thr / 2).max(1).min(sh.len());
    eprintln!("shards {} ({}) + unmapped, readers {nr}, htslib threads each {}", sh.len() - 1, if bed.is_some() { "bed regions" } else { "10 Mb windows" }, (thr / nr).max(1));
    let (next, ci, nrec, prx) = (AtomicUsize::new(0), AtomicU64::new(0), AtomicU64::new(0), Mutex::new(prx));
    std::thread::scope(|sc| for _ in 0..nr { sc.spawn(|| {
        let mut r = open(thr / nr); let mut rec = bam::Record::new(); let mut n = 0u64;
        let newbuf = || prx.lock().unwrap().try_recv().unwrap_or_else(|_| Vec::with_capacity(chunk() + (1 << 17)));
        let mut buf = newbuf(); buf.clear();
        loop {
            let i = next.fetch_add(1, Relaxed); if i >= sh.len() { break; }
            let (t, a, pe) = match sh[i] { Some((t, a, e, pe)) => { r.fetch((t as i32, a, e)).unwrap_or_else(|x| die(&format!("fetch {}:{a}-{e}: {x}", String::from_utf8_lossy(&names[t as usize])))); (t as i32, a, pe) }
                None => { r.fetch(bam::FetchDefinition::Unmapped).unwrap_or_else(|x| die(&format!("fetch unmapped: {x}"))); (-1, i64::MIN, i64::MIN) } };
            while let Some(x) = r.read(&mut rec) {
                x.unwrap_or_else(|e| die(&format!("read {input}: {e}")));
                let (fl, p) = (rec.flags(), rec.pos()); if fl & 0x900 != 0 || (t >= 0 && p < a && p < pe) { continue; }   // secondary/supplementary, or seen in the previous shard
                if buf.is_empty() { buf.push(0xFF); }   // binary record: flag u16, tid i32, pos i64, mapq u8, l_seq u32, 4-bit seq
                let l = rec.seq_len(); buf.extend_from_slice(&fl.to_le_bytes()); buf.extend_from_slice(&rec.tid().to_le_bytes()); buf.extend_from_slice(&p.to_le_bytes()); buf.push(rec.mapq()); buf.extend_from_slice(&(l as u32).to_le_bytes());
                buf.extend_from_slice(&rec.seq().encoded[..(l + 1) / 2]); n += 1;
                if buf.len() >= chunk() { let mut b = newbuf(); b.clear(); tx.send((ci.fetch_add(1, Relaxed), std::mem::replace(&mut buf, b))).unwrap(); }
            }
        }
        if buf.len() > 1 { tx.send((ci.fetch_add(1, Relaxed), buf)).unwrap(); }
        nrec.fetch_add(n, Relaxed);
    }); });
    let n = nrec.load(Relaxed); eprintln!("records read (primary) {n}"); n
}

fn count(src: Src, fa_path: &str, name: &str, bin: u64, thr: usize, out: &str) {
    let t0 = Instant::now();
    // reader: stdin read during the marker load (up to 64 chunks of 8 MB queued), chunks cut at line ends
    let q: usize = 32; eprintln!("queue {q} chunks of {} MB", chunk() >> 20);
    let (tx, rx) = mpsc::sync_channel::<(u64, Vec<u8>)>(q); let rx = Arc::new(Mutex::new(rx));
    let (ptx, prx) = mpsc::channel::<Vec<u8>>();   // used chunk buffers back to the reader (no new pages per chunk)
    let reader = spawn_reader(src, tx, prx, thr);
    let mut mk = load(fa_path, thr); let ncl = mk.names.len(); let n = mk.key.len();
    let per: Vec<usize> = (0..ncl).map(|j| mk.ocls.iter().filter(|&&c| c as usize == j).count()).collect();
    let tab = std::mem::replace(&mut mk.tab, Table { s: Vec::new(), cls: Vec::new(), mask: 0 });
    eprintln!("markers lines {} distinct {} classes {:?} per class {:?}; table {} MB; load {:.1} s; {}", mk.order.len(), n, mk.names, per, (tab.s.len() * 16) >> 20, t0.elapsed().as_secs_f64(), rss());
    let bpk: usize = 16;   // Bloom bits per key
    let bloom = { let mut b = Bloom::new(n, bpk.max(1)); for &x in &mk.key { b.insert(hk(x)); } Arc::new(b) };
    eprintln!("bloom {} MB ({} bits/key)", bloom.b.len() * 8 >> 20, bpk.max(1));
    let tab = Arc::new(tab);
    let cb: u32 = 15;   // per-thread cache: 2^cb entries
    let (ftx, frx) = mpsc::sync_channel::<(u64, Vec<u8>)>(thr * 4);
    std::fs::create_dir_all(out).unwrap_or_else(|e| die(&format!("{out}: {e}")));
    let any_path = format!("{out}/any.fa");
    let writer = std::thread::Builder::new().name("yr-write".into()).spawn(move || { let mut w = BufWriter::with_capacity(1 << 22, std::fs::File::create(&any_path).unwrap_or_else(|e| die(&format!("{any_path}: {e}")))); let mut next = 0u64; let mut held: std::collections::BTreeMap<u64, Vec<u8>> = std::collections::BTreeMap::new();   // marker reads in input order
        for (i, b) in frx { held.insert(i, b); while let Some(b) = held.remove(&next) { w.write_all(&b).unwrap_or_else(|e| die(&format!("write any.fa: {e}"))); next += 1; } }
        if !held.is_empty() { die("marker reads: chunks missing"); } w.flush().unwrap(); }).unwrap();
    let workers: Vec<_> = (0..thr).map(|_| {
        let (rx, ftx, ptx, tab, bloom) = (rx.clone(), ftx.clone(), ptx.clone(), tab.clone(), bloom.clone());
        std::thread::Builder::new().name("yr-count".into()).spawn(move || {
            let mut loc = Local { contigs: HashMap::new(), cname: vec![], bins: HashMap::new(), rk: vec![], rv: vec![], reads: 0, kmers: 0, hits: 0, cmiss: 0, bpos: 0, cache: vec![Ent { key: EMPTY, sc: 0, pend: 0, cls: 0 }; 1 << cb] }; let mut kb: Vec<(u64, u64, bool)> = Vec::with_capacity(1024); let mut sq: Vec<u8> = Vec::with_capacity(1 << 16); let mut cx = Ctx { bloom: &bloom, tab: &tab, ncl, bin, h: vec![0; ncl], last: (Vec::new(), u64::MAX, 0) };
            loop {
                let c = { let r = rx.lock().unwrap(); r.recv() };
                let Ok((ci, c)) = c else { break };
                let mut fa = Vec::with_capacity(1 << 16);
                if c.first() == Some(&0xFF) { work_bin(&c, &mut loc, &mut cx, &mut fa, &mut kb, &mut sq); } else { work(&c, &mut loc, &mut cx, &mut fa, &mut kb); }
                ftx.send((ci, fa)).unwrap(); let _ = ptx.send(c);   // fails only after the reader ended: buffer then freed
            }
            for e in &loc.cache { evict(e, &tab); }
            loc.bins = HashMap::new(); loc.cache = vec![]; loc   // lookup maps not needed after counting
        })
    }).collect();
    drop(ftx);
    let bytes_in = reader.join().unwrap();
    let locs: Vec<Local> = workers.into_iter().map(|w| w.unwrap().join().unwrap()).collect();
    writer.join().unwrap();
    let t_count = t0.elapsed().as_secs_f64(); eprintln!("counting done; {}", rss());
    let tp = Instant::now(); let lap = |w: &str| eprintln!("post {w} at {:.2} s", tp.elapsed().as_secs_f64());
    // bins: global contig ids, rows sorted by (contig name, bin); worker rows summed while writing (no merged copy)
    let w = 3 + 3 * ncl; let (mut gname, mut gid): (Vec<Vec<u8>>, HashMap<Vec<u8>, u32>) = (vec![], HashMap::new());
    let (mut reads, mut kmers, mut hits, mut cmiss, mut bpos) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut ord: Vec<(u32, u64, u32, u32)> = Vec::with_capacity(locs.iter().map(|l| l.rk.len()).sum());   // (contig, bin, worker, row)
    for (wi, l) in locs.iter().enumerate() {
        reads += l.reads; kmers += l.kmers; cmiss += l.cmiss; bpos += l.bpos; hits += l.hits;
        let g: Vec<u32> = l.cname.iter().map(|c| *gid.entry(c.clone()).or_insert_with(|| { gname.push(c.clone()); gname.len() as u32 - 1 })).collect();
        for (r, &(ci, b)) in l.rk.iter().enumerate() { ord.push((g[ci as usize], b, wi as u32, r as u32)); }
    }
    if reads == 0 { die("0 reads"); }
    let mut rank = vec![0u32; gname.len()]; { let mut o: Vec<usize> = (0..gname.len()).collect(); o.sort_by(|&a, &b| gname[a].cmp(&gname[b])); for (r, &i) in o.iter().enumerate() { rank[i] = r as u32; } }
    ord.sort_unstable_by_key(|&(c, b, wi, r)| (rank[c as usize], b, wi, r));
    let nb = ord.len() - ord.windows(2).filter(|p| p[0].0 == p[1].0 && p[0].1 == p[1].1).count();
    eprintln!("bin records over workers {} (merged {nb}); {}", ord.len(), rss());
    eprintln!("reads {reads} bins {nb} bytes {} k-mers {kmers} cache misses {cmiss} ({:.2}%) bloom positives {bpos} hits {hits}; counting {:.1} s", bytes_in, 100.0 * cmiss as f64 / kmers.max(1) as f64, t_count);
    lap("bins merged");
    let mut cnt: Vec<u64> = vec![0; n]; for sl in tab.s.iter() { if sl.key != EMPTY { cnt[sl.idx as usize] = sl.cnt.load(Relaxed) as u64; } }
    drop(tab); drop(bloom);   // counts extracted: table and filter freed
    // per class lists in marker order
    let mut cv: Vec<Vec<u64>> = vec![vec![]; ncl]; let mut gv: Vec<Vec<u8>> = vec![vec![]; ncl];
    for (k, &i) in mk.order.iter().enumerate() { let c = mk.ocls[k] as usize; cv[c].push(cnt[i as usize]); gv[c].push(mk.gc[k]); }
    // check: bins hits = counts per class
    for j in 0..ncl { let s: u64 = cv[j].iter().sum(); let bh: u64 = locs.iter().map(|l| l.rv.iter().skip(3 + 3 * j).step_by(w).sum::<u64>()).sum(); eprintln!("check {} counts {s} bins {bh}", mk.names[j]); if s != bh { die(&format!("{} bins hits {bh} vs counts {s}", mk.names[j])); } }
    let wr = |p: &str| BufWriter::with_capacity(1 << 22, std::fs::File::create(format!("{out}/{p}")).unwrap_or_else(|e| die(&format!("{p}: {e}"))));
    let t1 = Instant::now();
    { // counts.tsv: rounds of thr parallel blocks (1M lines each), written in order
        let (key, cnt, per) = (&mk.key, &cnt, 1usize << 20); let mut o = wr("counts.tsv");
        for r0 in (0..n).step_by(per * thr) {
            let parts: Vec<Vec<u8>> = std::thread::scope(|sc| { let hs: Vec<_> = (r0..(r0 + per * thr).min(n)).step_by(per).map(|a| sc.spawn(move || { let mut b = Vec::with_capacity(per * 40); for i in a..(a + per).min(n) { b.extend_from_slice(&decode(key[i])); b.push(b'\t'); itoa(&mut b, cnt[i]); b.push(b'\n'); } b })).collect(); hs.into_iter().map(|h| h.join().unwrap()).collect() });
            for b in &parts { o.write_all(b).unwrap_or_else(|e| die(&format!("write counts.tsv: {e}"))); }
        }
        o.flush().unwrap();
    }
    eprintln!("write counts {:.2} s", t1.elapsed().as_secs_f64());
    {
        let mut o = wr("bins.tsv"); let mut v = vec![0u64; w]; let mut b: Vec<u8> = Vec::with_capacity(1 << 20);
        write!(o, "contig\tstart\treads\treads_mq20\tduplicates").unwrap();
        for c in &mk.names { write!(o, "\t{c}_hits\t{c}_hits_mq20\t{c}_reads").unwrap(); } writeln!(o).unwrap();
        for (i, &(c, bn, wi, r)) in ord.iter().enumerate() {
            let x = &locs[wi as usize].rv[r as usize * w..(r as usize + 1) * w]; for (a, y) in v.iter_mut().zip(x) { *a += y; }
            if ord.get(i + 1).map_or(true, |n| n.0 != c || n.1 != bn) {
                b.extend_from_slice(&gname[c as usize]); b.push(b'\t'); itoa(&mut b, bn * bin); for a in v.iter_mut() { b.push(b'\t'); itoa(&mut b, *a); *a = 0; } b.push(b'\n');
                if b.len() >= 1 << 20 { o.write_all(&b).unwrap_or_else(|e| die(&format!("write bins.tsv: {e}"))); b.clear(); }
            }
        }
        o.write_all(&b).unwrap_or_else(|e| die(&format!("write bins.tsv: {e}"))); o.flush().unwrap();
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
    eprintln!("y_dose {:?} ysinglecopy_depth {:?} chr2depth_mean {:?}; total {:.1} s; {}", 2.0 * xd / ad, xd, ad, t0.elapsed().as_secs_f64(), rss());
}

const TECH: [&str; 5] = ["illumina", "element", "hifi", "ont", "hic"];
// estimate: length = K * sum_gc(class) / sum(ysinglecopy), K of technology x (model column 3, empty = illumina)
fn estimate(summary: &str, model: &str, x: &str, o: &mut dyn Write) {
    let rd = |p: &str| std::fs::read_to_string(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let (sm, md) = (rd(summary), rd(model));
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
    for l in md.lines().filter(|l| !l.starts_with('#')) { let f: Vec<&str> = l.split('\t').collect(); if f.len() >= 2 && f[0] != "class" && f.get(2).map_or("illumina", |t| t.trim()) == x { k.insert(f[0].to_string(), f[1].parse().unwrap_or_else(|_| die(&format!("model K {}", f[1])))); } }
    if k.is_empty() { die(&format!("model {model} has no K for technology {x}")); }
    if x == "element" || x == "hic" { eprintln!("WARNING technology {x} is experimental"); }
    writeln!(o, "sample\tclass\tlength_bp\tmarker_sum_gc\tysinglecopy_sum\ty_dose").unwrap();
    let mut n = 0;
    for (s, c, v) in &sg { if let Some(kk) = k.get(c) { writeln!(o, "{s}\t{c}\t{:.0}\t{v:.1}\t{ysum:.0}\t{dose}", kk * v / ysum).unwrap(); n += 1; } }
    for c in k.keys() { if !sg.iter().any(|t| &t.1 == c) { die(&format!("model class {c} missing in summary")); } }
    eprintln!("estimate technology {x} classes {n} y_dose {dose}");
}

const USAGE: &str = "usage:
  yruler -i <cram|bam|-> -m <markers> -o <outdir> [-r <ref.fa>] [-b <regions.bed>] [-M <model.tsv>] [-x <tech>] [-n <name>] [-t <threads>]
  yruler estimate <summary.tsv> <model.tsv> [-x <tech>]
-i -: SAM on stdin
-x: illumina (default), element (experimental), hifi, ont, hic (experimental)";
fn tech(t: &str) -> &'static str { TECH.iter().find(|&&x| x == t).copied().unwrap_or_else(|| die(&format!("-x {t}: use one of {}", TECH.join(", ")))) }
fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    // static libcurl: point it to the system CA bundle
    if std::env::var_os("CURL_CA_BUNDLE").is_none() { if let Some(p) = ["/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt", "/etc/ssl/cert.pem"].iter().find(|p| std::path::Path::new(p).exists()) { unsafe { std::env::set_var("CURL_CA_BUNDLE", p) }; } }
    match a.first().map(|s| s.as_str()) {
        Some("estimate") => { let x = match a.len() { 3 => "illumina", 5 if a[3] == "-x" => tech(&a[4]), _ => die(USAGE) }; estimate(&a[1], &a[2], x, &mut io::stdout().lock()); return; }
        Some("--version") | Some("-V") => { println!("yruler {}", env!("CARGO_PKG_VERSION")); return; }
        None | Some("-h") | Some("--help") => { eprintln!("{USAGE}"); std::process::exit(if a.is_empty() { 1 } else { 0 }); }
        _ => {}
    }
    let mut o: HashMap<&str, String> = HashMap::new(); let mut i = 0;
    while i < a.len() { let k = a[i].as_str(); if !matches!(k, "-i" | "-m" | "-o" | "-r" | "-b" | "-M" | "-x" | "-n" | "-t" | "--bin") { die(&format!("unknown option {k}\n{USAGE}")); } let v = a.get(i + 1).unwrap_or_else(|| die(&format!("{k} needs a value"))); o.insert(k, v.clone()); i += 2; }
    let need = |k: &str| o.get(k).cloned().unwrap_or_else(|| die(&format!("missing {k}\n{USAGE}")));
    let (input, markers, out) = (need("-i"), need("-m"), need("-o")); let x = o.get("-x").map_or("illumina", |t| tech(t));
    let thr = o.get("-t").map_or(std::thread::available_parallelism().map_or(1, |n| n.get()), |t| t.parse().unwrap_or_else(|_| die("-t")));
    let bin: u64 = o.get("--bin").map_or(500, |b| b.parse().unwrap_or_else(|_| die("--bin")));
    let name = o.get("-n").cloned().unwrap_or_else(|| { let b = std::path::Path::new(&input).file_name().map_or("sample".into(), |f| f.to_string_lossy().to_string()); b.split('.').next().unwrap_or("sample").to_string() });
    let src = if input == "-" { Src::Stdin } else { Src::Hts { input: input.clone(), reference: o.get("-r").cloned(), bed: o.get("-b").cloned() } };
    eprintln!("yruler {} input {input} markers {markers} name {name} threads {thr} technology {x} regions {}", env!("CARGO_PKG_VERSION"), o.get("-b").map_or("whole file", |s| s.as_str()));
    count(src, &markers, &name, bin, thr, &out);
    if let Some(m) = o.get("-M") { let p = format!("{out}/lengths.tsv"); let mut w = BufWriter::new(std::fs::File::create(&p).unwrap_or_else(|e| die(&format!("{p}: {e}")))); estimate(&format!("{out}/summary.tsv"), m, x, &mut w); w.flush().unwrap(); print!("{}", std::fs::read_to_string(&p).unwrap()); }
}

#[cfg(test)]
mod tests {
    #[test]
    fn masked_superstring() {   // uppercase starts a k-mer; header name before '_'
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGTa";
        let mut ms = b">c1_x\n".to_vec(); ms.extend(s.iter().enumerate().map(|(i, b)| if i < 2 { b.to_ascii_uppercase() } else { b.to_ascii_lowercase() })); ms.push(b'\n');
        let r = super::ms_records(&ms); assert_eq!(r.len(), 1); assert_eq!(r[0].0, b"c1");
        let st: Vec<usize> = (0..r[0].1.len() + 1 - super::K).filter(|&i| r[0].1[i].is_ascii_uppercase()).collect(); assert_eq!(st, vec![0, 1]);
        assert_eq!(super::canon(&r[0].1[1..32]), super::canon(b"CGTACGTACGTACGTACGTACGTACGTACGT"));
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

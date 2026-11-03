// yruler: Y satellite array lengths from reads
// count: SAM on stdin -> marker k-mer counts per class
// estimate: class sums -> array lengths, Y dosage
use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering::Relaxed};
use std::sync::{mpsc, Arc, Mutex};
use std::time::Instant;

const K: usize = 31;
const MASK: u64 = (1u64 << 62) - 1;
const EMPTY: u64 = u64::MAX;
const CHUNK: usize = 2 << 20; // read chunk bytes (2 MB: same speed as 8 to 64, least RSS); a longer line grows its chunk
fn chunk() -> usize {
    CHUNK
}

static LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
};

#[inline(always)]
fn kmer_hash(x: u64) -> u64 {
    let m = x.wrapping_mul(0x9e3779b97f4a7c15);
    m ^ (m >> 29)
} // k-mer hash: one multiply, high bits folded down

fn die(m: &str) -> ! {
    eprintln!("ERROR yruler: {m}");
    std::process::exit(1)
}
// resident and peak memory from /proc (Linux), for the stage log
fn rss() -> String {
    std::fs::read_to_string("/proc/self/status").map_or("rss n/a".into(), |t| {
        t.lines()
            .filter(|l| l.starts_with("VmRSS") || l.starts_with("VmHWM"))
            .map(|l| l.split_whitespace().take(2).collect::<Vec<_>>().join(" "))
            .collect::<Vec<_>>()
            .join(" ")
    })
}

// word Bloom filter, checked only on k-mer cache misses: 4 bits in one 64-bit word per key (one load, one AND)
struct Bloom {
    b: Vec<u64>,
    mask: u64,
}
impl Bloom {
    fn new(n: usize, bits_per_key: usize) -> Self {
        let nw = ((n * bits_per_key + 63) / 64).next_power_of_two();
        let mut b: Vec<u64> = Vec::with_capacity(nw);
        hugepage(&b);
        b.resize(nw, 0);
        Bloom { b, mask: nw as u64 - 1 }
    }
    #[inline(always)]
    fn pat(h: u64) -> u64 {
        (1u64 << (h & 63)) | (1u64 << ((h >> 6) & 63)) | (1u64 << ((h >> 12) & 63)) | (1u64 << ((h >> 18) & 63))
    }
    #[inline(always)]
    fn at(&self, h: u64) -> *const u64 {
        unsafe { self.b.as_ptr().add(((h >> 32) & self.mask) as usize) }
    }
    fn insert(&mut self, h: u64) {
        let i = ((h >> 32) & self.mask) as usize;
        self.b[i] |= Self::pat(h);
    }
    #[inline(always)]
    fn has(&self, h: u64) -> bool {
        let m = Self::pat(h);
        unsafe { *self.at(h) & m == m }
    }
}

// open-addressing table: 16-byte slot (key, count, marker index), class bits per slot (up to 16 classes)
#[repr(C, align(16))]
struct Slot {
    key: u64,
    cnt: AtomicU32,
    idx: u32,
}
fn hugepage<T>(v: &Vec<T>) {
    // ask the kernel for transparent huge pages before the first touch (fewer TLB misses)
    let (p, n) = (v.as_ptr() as usize, v.capacity() * std::mem::size_of::<T>());
    let a = (p + (1 << 21) - 1) & !((1 << 21) - 1);
    if n > (4 << 21) {
        unsafe {
            libc::madvise(
                a as *mut libc::c_void,
                (p + n - a) & !((1 << 21) - 1),
                libc::MADV_HUGEPAGE,
            );
        }
    }
}
struct Table {
    s: Vec<Slot>,
    cls: Vec<u16>,
    mask: usize,
}
impl Table {
    fn new(n: usize, thr: usize) -> Self {
        // slots written in parallel (first touch spread over threads)
        let s = (n * 2).next_power_of_two();
        let mut v: Vec<Slot> = Vec::with_capacity(s);
        hugepage(&v);
        let (p, per) = (v.as_mut_ptr() as usize, s.div_ceil(thr));
        std::thread::scope(|sc| {
            for t in 0..thr {
                sc.spawn(move || {
                    let p = p as *mut Slot;
                    for i in t * per..((t + 1) * per).min(s) {
                        unsafe {
                            p.add(i).write(Slot {
                                key: EMPTY,
                                cnt: AtomicU32::new(0),
                                idx: 0,
                            });
                        }
                    }
                });
            }
        });
        unsafe {
            v.set_len(s);
        }
        Table {
            s: v,
            cls: vec![0; s],
            mask: s - 1,
        }
    }
    // returns the marker index of key, inserting it with index val; class bit c is ORed in
    fn insert(&mut self, key: u64, h: u64, val: u32, c: usize) -> (u32, bool) {
        let mut i = (h >> 7) as usize & self.mask;
        loop {
            let sl = &mut self.s[i];
            if sl.key == EMPTY {
                sl.key = key;
                sl.idx = val;
                self.cls[i] |= 1 << c;
                return (val, true);
            }
            if sl.key == key {
                self.cls[i] |= 1 << c;
                return (sl.idx, false);
            }
            i = (i + 1) & self.mask;
        }
    }
    #[inline(always)]
    fn at(&self, h: u64) -> *const Slot {
        unsafe { self.s.as_ptr().add((h >> 7) as usize & self.mask) }
    }
    #[inline(always)]
    fn find(&self, key: u64, h: u64) -> Option<usize> {
        let mut i = (h >> 7) as usize & self.mask;
        loop {
            let sl = &self.s[i];
            if sl.key == key {
                return Some(i);
            }
            if sl.key == EMPTY {
                return None;
            }
            i = (i + 1) & self.mask;
        }
    }
}
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch(p as *const i8, std::arch::x86_64::_MM_HINT_T0);
    }
}

fn canon(s: &[u8]) -> Option<u64> {
    if s.len() != K {
        return None;
    }
    let (mut f, mut r) = (0u64, 0u64);
    for &b in s {
        let c = LUT[b as usize];
        if c > 3 {
            return None;
        }
        f = (f << 2) | c as u64;
        r = (r >> 2) | ((3 - c as u64) << 60);
    }
    Some(f.min(r))
}
fn decode(x: u64) -> [u8; K] {
    let mut o = [0u8; K];
    for i in 0..K {
        o[i] = b"ACGT"[((x >> (2 * (K - 1 - i))) & 3) as usize];
    }
    o
}

struct Markers {
    key: Vec<u64>,
    order: Vec<u32>,
    ocls: Vec<u8>,
    gc: Vec<u8>,
    names: Vec<String>,
    tab: Table,
}
// class name of a header line: text before the first '_' (trailing whitespace trimmed first)
fn class_name(h: &[u8]) -> &[u8] {
    let mut e = h.len();
    while e > 0 && h[e - 1].is_ascii_whitespace() {
        e -= 1;
    }
    let h = &h[..e];
    &h[..h.iter().position(|&b| b == b'_').unwrap_or(h.len())]
}
// masked superstring records (class name, sequence): uppercase base = a k-mer starts there
fn ms_records(raw: &[u8]) -> Vec<(&[u8], Vec<u8>)> {
    let mut r: Vec<(&[u8], Vec<u8>)> = vec![];
    for l in raw.split(|&b| b == b'\n') {
        let l = l.strip_suffix(b"\r").unwrap_or(l);
        if l.first() == Some(&b'>') {
            r.push((class_name(&l[1..]), vec![]));
        } else if let Some(x) = r.last_mut() {
            x.1.extend_from_slice(l);
        } else if !l.is_empty() {
            die("sequence before the first header");
        }
    }
    r
}
type Part<'a> = (Vec<(u64, u8, u8)>, Vec<&'a [u8]>); // (key, G+C, local header index or 255), header names
fn gc_count(w: &[u8]) -> u8 {
    w.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count() as u8
}
// k-mers in file order, parsed in parallel pieces; first occurrence = marker index
fn load(path: &str, thr: usize) -> Markers {
    let raw = std::fs::read(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    let raw = if path.ends_with(".zst") {
        zstd::stream::decode_all(&raw[..]).unwrap_or_else(|e| die(&format!("{path}: zstd {e}")))
    } else {
        raw
    };
    // one k-mer per line, else masked superstring
    let ms = raw
        .split(|&b| b == b'\n')
        .find(|l| !l.is_empty() && l[0] != b'>')
        .map_or(false, |l| {
            let l = l.strip_suffix(b"\r").unwrap_or(l);
            l.len() != K || l.iter().any(|b| b.is_ascii_lowercase())
        });
    let recs = if ms { ms_records(&raw) } else { vec![] };
    let parts: Vec<Part> = if ms {
        // record pieces of about starts/thr k-mer starts
        let ns: usize = recs.iter().map(|r| (r.1.len() + 1).saturating_sub(K)).sum();
        let per = (ns / thr).max(1 << 16);
        let tasks: Vec<(usize, usize, usize)> = recs
            .iter()
            .enumerate()
            .flat_map(|(ri, r)| {
                let n = (r.1.len() + 1).saturating_sub(K);
                (0..n).step_by(per).map(move |a| (ri, a, (a + per).min(n)))
            })
            .collect();
        let p: Vec<Part> = std::thread::scope(|sc| {
            let recs = &recs;
            let hs: Vec<_> = tasks
                .iter()
                .map(|&(ri, a, b)| {
                    sc.spawn(move || {
                        let s = &recs[ri].1;
                        let mut v = Vec::with_capacity((b - a) / 2);
                        for i in a..b {
                            if s[i].is_ascii_uppercase() {
                                let w = &s[i..i + K];
                                v.push((
                                    canon(w).unwrap_or_else(|| die("masked superstring k-mer with non-ACGT")),
                                    gc_count(w),
                                    0,
                                ));
                            }
                        }
                        (v, vec![recs[ri].0])
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        eprintln!(
            "markers: masked superstring {} bytes, {} records, {} k-mer starts, {} pieces",
            raw.len(),
            recs.len(),
            ns,
            tasks.len()
        );
        p
    } else {
        // file cut at line ends into thr pieces
        let n = raw.len();
        let mut cuts = vec![0usize];
        for t in 1..thr {
            let mut p = (n * t / thr).max(*cuts.last().unwrap());
            while p < n && p > 0 && raw[p - 1] != b'\n' {
                p += 1;
            }
            cuts.push(p);
        }
        cuts.push(n);
        std::thread::scope(|sc| {
            let hs: Vec<_> = cuts
                .windows(2)
                .map(|w| {
                    let raw = &raw[w[0]..w[1]];
                    sc.spawn(move || {
                        let (mut v, mut nm, mut cur): (Vec<(u64, u8, u8)>, Vec<&[u8]>, u8) =
                            (Vec::with_capacity(raw.len() / 64 + 16), vec![], 255);
                        for l in raw.split(|&b| b == b'\n') {
                            if l.is_empty() {
                                continue;
                            }
                            if l[0] == b'>' {
                                let c = class_name(&l[1..]);
                                if nm.last() != Some(&c) {
                                    if nm.len() == 254 {
                                        die("more than 254 header names in one piece");
                                    }
                                    nm.push(c);
                                }
                                cur = (nm.len() - 1) as u8;
                                continue;
                            }
                            let l = l.strip_suffix(b"\r").unwrap_or(l);
                            v.push((
                                canon(l)
                                    .unwrap_or_else(|| die(&format!("marker line of length {} or non-ACGT", l.len()))),
                                gc_count(l),
                                cur,
                            ));
                        }
                        (v, nm)
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        })
    };
    let mut names: Vec<String> = vec![];
    let mut cur = usize::MAX;
    let tot: usize = parts.iter().map(|p| p.0.len()).sum();
    if tot >= u32::MAX as usize {
        die("more than 2^32 - 1 markers");
    }
    let mut tab = Table::new(tot, thr);
    let (mut key, mut order, mut ocls, mut gc) = (
        vec![],
        Vec::with_capacity(tot),
        Vec::with_capacity(tot),
        Vec::with_capacity(tot),
    );
    for (v, nm) in &parts {
        let gl: Vec<usize> = nm
            .iter()
            .map(|c| {
                let c = String::from_utf8_lossy(c).to_string();
                names.iter().position(|x| *x == c).unwrap_or_else(|| {
                    names.push(c);
                    names.len() - 1
                })
            })
            .collect();
        if names.len() > 16 {
            die("more than 16 classes");
        }
        for (k, &(x, g, lc)) in v.iter().enumerate() {
            if lc != 255 {
                cur = gl[lc as usize];
            } else if cur == usize::MAX {
                die("sequence before the first header");
            }
            if let Some(nx) = v.get(k + 16) {
                prefetch(tab.at(kmer_hash(nx.0)));
            }
            let (i, new) = tab.insert(x, kmer_hash(x), key.len() as u32, cur);
            if new {
                key.push(x);
            }
            order.push(i);
            ocls.push(cur as u8);
            gc.push(g);
        }
        if let Some(&c) = gl.last() {
            cur = c;
        }
    }
    Markers {
        key,
        order,
        ocls,
        gc,
        names,
        tab,
    }
}

// per-thread bin record (rk key, rv flat values): reads, reads MAPQ >= 20, duplicates, then per class hits, hits MAPQ >= 20, reads with a hit
struct Local {
    contigs: HashMap<Vec<u8>, u32>,
    cname: Vec<Vec<u8>>,
    bins: HashMap<(u32, u64), usize>,
    rk: Vec<(u32, u64)>,
    rv: Vec<u64>,
    reads: u64,
    kmers: u64,
    hits: u64,
    cmiss: u64,
    bpos: u64,
    cache: Vec<Ent>,
}
// per-thread cache of recent k-mers (sorted input: a k-mer comes back in the next reads); counts added to the table on eviction
#[derive(Clone, Copy)]
struct Ent {
    key: u64,
    sc: u32,
    pend: u16,
    cls: u16,
} // 16 bytes; sc = table slot + 1, 0 = not a marker
fn evict(e: &Ent, tab: &Table) {
    if e.sc != 0 && e.pend > 0 {
        tab.s[e.sc as usize - 1].cnt.fetch_add(e.pend as u32, Relaxed);
    }
}
#[inline(always)]
fn hit(e: &mut Ent, tab: &Table, h: &mut [u64]) {
    // count, flush before u16 overflow
    if e.pend == u16::MAX {
        tab.s[e.sc as usize - 1].cnt.fetch_add(e.pend as u32, Relaxed);
        e.pend = 0;
    }
    e.pend += 1;
    let mut m = e.cls;
    while m != 0 {
        h[m.trailing_zeros() as usize] += 1;
        m &= m - 1;
    }
}

#[inline(always)]
fn memchr(s: &[u8], c: u8) -> Option<usize> {
    let p = unsafe { libc::memchr(s.as_ptr() as *const libc::c_void, c as i32, s.len()) };
    if p.is_null() {
        None
    } else {
        Some(p as usize - s.as_ptr() as usize)
    }
} // glibc SIMD memchr
#[inline(always)]
fn field<'a>(l: &'a [u8], q: &mut usize) -> &'a [u8] {
    // next tab-separated field from q (last field: to line end)
    if *q > l.len() {
        die("SAM line with too few fields");
    }
    let r = &l[*q..];
    let e = memchr(r, b'\t').unwrap_or(r.len());
    *q += e + 1;
    &r[..e]
}
fn parse_u64(s: &[u8]) -> u64 {
    let mut v = 0u64;
    for &b in s {
        if !b.is_ascii_digit() {
            die("non-numeric SAM field");
        }
        v = v * 10 + (b - b'0') as u64;
    }
    v
}

fn itoa(b: &mut Vec<u8>, mut v: u64) {
    let mut t = [0u8; 20];
    let mut i = 20;
    loop {
        i -= 1;
        t[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    b.extend_from_slice(&t[i..]);
}
struct Ctx<'a> {
    bloom: &'a Bloom,
    tab: &'a Table,
    ncl: usize,
    bin: u64,
    h: Vec<u64>,
    last: (Vec<u8>, u64, usize),
}
// one read: k-mer hits, bin record, marker read
fn count_read(
    c: &mut Ctx,
    loc: &mut Local,
    fa: &mut Vec<u8>,
    kb: &mut Vec<(u64, u64, bool)>,
    rname: &[u8],
    pos: u64,
    mapq: u64,
    fl: u64,
    seq: &[u8],
    qual: &[u8],
) {
    let (w, b, mq) = (3 + 3 * c.ncl, pos.saturating_sub(1) / c.bin, mapq >= 20); // pos is 1-based
    let ri = if c.last.1 == b && c.last.0 == rname {
        c.last.2
    } else {
        let ci = match loc.contigs.get(rname) {
            Some(&x) => x,
            None => {
                let x = loc.cname.len() as u32;
                loc.contigs.insert(rname.to_vec(), x);
                loc.cname.push(rname.to_vec());
                x
            }
        };
        let ri = *loc.bins.entry((ci, b)).or_insert_with(|| {
            loc.rk.push((ci, b));
            loc.rv.resize(loc.rv.len() + w, 0);
            loc.rk.len() - 1
        });
        c.last = (rname.to_vec(), b, ri);
        ri
    };
    loc.reads += 1;
    c.h.iter_mut().for_each(|x| *x = 0);
    let h = &mut c.h;
    if seq.contains(&b'=') {
        EQB.fetch_add(seq.iter().filter(|&&b| b == b'=').count() as u64, Relaxed);
    }
    let (mut f, mut r, mut n, mut n_k) = (0u64, 0u64, 0usize, 0u64);
    kb.clear();
    let cm = loc.cache.len() - 1;
    let qm = QUAL.get().filter(|_| qual.len() == seq.len());
    let seed = kmer_hash(fl ^ (pos << 16) ^ ((seq.len() as u64) << 48));
    let (mut s, mut bad) = (0i64, usize::MAX);
    if seq != b"*" {
        for (i, &ch) in seq.iter().enumerate() {
            let ch = LUT[ch as usize];
            if ch > 3 {
                n = 0;
                f = 0;
                r = 0;
                s = 0;
                continue;
            }
            f = ((f << 2) | ch as u64) & MASK;
            r = (r >> 2) | ((3 - ch as u64) << 60);
            n += 1;
            if let Some((qmin, lp)) = qm {
                if *qmin > 0 {
                    if qual[i] < *qmin {
                        bad = i;
                    }
                } else {
                    s += lp[qual[i] as usize];
                    if n > K {
                        s -= lp[qual[i - K] as usize];
                    }
                }
            }
            if n >= K {
                n_k += 1;
                if let Some((qmin, _)) = qm {
                    // low-quality k-mer: dropped (min Q) or kept with P(all 31 bases right)
                    if *qmin > 0 {
                        if bad != usize::MAX && bad + K > i {
                            continue;
                        }
                    } else if (splitmix(seed ^ f ^ ((i as u64) << 40)) >> 11) as f64 / (1u64 << 53) as f64
                        >= (s as f64 * 1e-9).exp()
                    {
                        continue;
                    }
                }
                let x = f.min(r);
                let hx = kmer_hash(x);
                let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
                if e.key == x {
                    if e.sc != 0 {
                        loc.hits += 1;
                        hit(e, c.tab, h);
                    }
                } else {
                    prefetch(c.bloom.at(hx));
                    kb.push((x, hx, false));
                }
            }
        }
    }
    loc.kmers += n_k;
    loc.cmiss += kb.len() as u64;
    for k in kb.iter_mut() {
        if c.bloom.has(k.1) {
            k.2 = true;
            prefetch(c.tab.at(k.1));
        }
    } // Bloom positives: prefetch slot
    for &(x, hx, bp) in kb.iter() {
        // cache misses: evict, look up, fill
        let e = unsafe { loc.cache.get_unchecked_mut((hx >> 40) as usize & cm) };
        if e.key != x {
            evict(e, c.tab);
            *e = match if bp {
                loc.bpos += 1;
                c.tab.find(x, hx)
            } else {
                None
            } {
                Some(i) => Ent {
                    key: x,
                    sc: i as u32 + 1,
                    pend: 0,
                    cls: c.tab.cls[i],
                },
                None => Ent {
                    key: x,
                    sc: 0,
                    pend: 0,
                    cls: 0,
                },
            };
        }
        if e.sc != 0 {
            loc.hits += 1;
            hit(e, c.tab, h);
        }
    }
    let rec = &mut loc.rv[ri * w..(ri + 1) * w];
    rec[0] += 1;
    if mq {
        rec[1] += 1;
    }
    if (fl >> 10) & 1 == 1 {
        rec[2] += 1;
    }
    let mut any = false;
    for j in 0..c.ncl {
        if h[j] > 0 {
            any = true;
            rec[3 + 3 * j] += h[j];
            if mq {
                rec[4 + 3 * j] += h[j];
            }
            rec[5 + 3 * j] += 1;
        }
    }
    if any && READS.get().is_some() {
        fa.push(b'>');
        fa.extend_from_slice(rname);
        fa.push(b'|');
        itoa(fa, b);
        fa.push(b'|');
        itoa(fa, mapq);
        fa.push(b'|');
        itoa(fa, fl);
        fa.push(b'\n');
        fa.extend_from_slice(seq);
        fa.push(b'\n');
    }
}
// SAM text chunk
fn work(chunk: &[u8], loc: &mut Local, c: &mut Ctx, fa: &mut Vec<u8>, kb: &mut Vec<(u64, u64, bool)>) {
    let mut p0 = 0;
    while p0 < chunk.len() {
        let e0 = memchr(&chunk[p0..], b'\n').map_or(chunk.len(), |i| p0 + i);
        let line = &chunk[p0..e0];
        p0 = e0 + 1;
        if line.is_empty() || line[0] == b'@' {
            continue;
        }
        let mut q = 0;
        let it = &mut q;
        field(line, it);
        let fl = parse_u64(field(line, it));
        let rname = field(line, it);
        let pos = parse_u64(field(line, it));
        let mapq = parse_u64(field(line, it));
        if fl & 0x900 != 0 {
            continue;
        } // secondary/supplementary, as for BAM/CRAM
        for _ in 0..4 {
            field(line, it);
        }
        let seq = field(line, it);
        let qv: Vec<u8> = if QUAL.get().is_some() {
            let q = field(line, it);
            if q == b"*" {
                vec![]
            } else {
                q.iter().map(|&b| b.wrapping_sub(33)).collect()
            }
        } else {
            vec![]
        };
        count_read(c, loc, fa, kb, rname, pos, mapq, fl, seq, &qv);
    }
}
// binary chunk (first byte 0xFF): flag u16, tid i32, pos i64, mapq u8, l_seq u32, 4-bit seq
fn work_bin(
    ch: &[u8],
    loc: &mut Local,
    c: &mut Ctx,
    fa: &mut Vec<u8>,
    kb: &mut Vec<(u64, u64, bool)>,
    sq: &mut Vec<u8>,
) {
    const NT: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let names = NAMES.get().unwrap_or_else(|| die("contig names not set"));
    let mut p = 1;
    while p < ch.len() {
        let fl = u16::from_le_bytes([ch[p], ch[p + 1]]) as u64;
        let tid = i32::from_le_bytes(ch[p + 2..p + 6].try_into().unwrap());
        let pos = (i64::from_le_bytes(ch[p + 6..p + 14].try_into().unwrap()) + 1).max(0) as u64;
        let mq = ch[p + 14] as u64;
        let n = u32::from_le_bytes(ch[p + 15..p + 19].try_into().unwrap()) as usize;
        p += 19;
        let e = &ch[p..p + (n + 1) / 2];
        p += (n + 1) / 2;
        sq.clear();
        if n == 0 {
            sq.push(b'*');
        } else {
            sq.extend((0..n).map(|i| NT[((e[i >> 1] >> (4 * (1 - (i & 1)))) & 15) as usize]));
        }
        let qs = if QUAL.get().is_some() {
            let q = &ch[p..p + n];
            p += n;
            q
        } else {
            &ch[p..p]
        };
        count_read(
            c,
            loc,
            fa,
            kb,
            if tid < 0 { b"*" } else { &names[tid as usize] },
            pos,
            mq,
            fl,
            sq,
            qs,
        );
    }
}

// Python repr of a float: shortest digits, fixed for exponents -4..15, else d.ddde-XX
fn py_float(x: f64) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0.0".into()
        } else {
            "0.0".into()
        };
    }
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let e = format!("{x:e}");
    let (m, ex) = e.split_once('e').unwrap();
    let ex: i32 = ex.parse().unwrap();
    let neg = m.starts_with('-');
    let d: String = m.trim_start_matches('-').chars().filter(|c| *c != '.').collect();
    let n = d.len() as i32;
    let s = if (-4..16).contains(&ex) {
        if ex >= 0 {
            if n <= ex + 1 {
                format!("{d}{}.0", "0".repeat((ex + 1 - n) as usize))
            } else {
                format!("{}.{}", &d[..(ex + 1) as usize], &d[(ex + 1) as usize..])
            }
        } else {
            format!("0.{}{d}", "0".repeat((-ex - 1) as usize))
        }
    } else {
        let mt = if n == 1 {
            d.clone()
        } else {
            format!("{}.{}", &d[..1], &d[1..])
        };
        format!("{mt}e{}{:02}", if ex < 0 { '-' } else { '+' }, ex.abs())
    };
    if neg {
        format!("-{s}")
    } else {
        s
    }
}

// GC curve from ysinglecopy, applied to array classes
fn median(v: &[u64]) -> (f64, String) {
    // statistics.median and its Python repr
    let mut s = v.to_vec();
    let n = s.len();
    let hi = *s.select_nth_unstable(n / 2).1; // O(n) selection, same values as a sort
    if n % 2 == 1 {
        (hi as f64, hi.to_string())
    } else {
        let lo = *s[..n / 2].iter().max().unwrap();
        let m = (lo + hi) as f64 / 2.0;
        (m, py_float(m))
    }
}
// G+C curve: for each G+C bin (number of G or C in the 31-mer), the mean ruler count in the bin / the mean ruler count.
// Counts above 10 x (mean + 1) are dropped. It needs 1e6 ruler counts, and 1,000 k-mers and 1e4 counts in a bin.
fn gc_bins(ruler: &[u64], ruler_by_gc: &[Vec<u64>]) -> Vec<(usize, f64)> {
    let total: u64 = ruler.iter().sum();
    let mean = total as f64 / ruler.len() as f64;
    let cap = 10.0 * (mean + 1.0);
    if (total as f64) < 1e6 || !(mean > 0.0) {
        return vec![];
    }
    let kept = |counts: &[u64]| {
        counts
            .iter()
            .filter(|&&x| x as f64 <= cap)
            .fold((0u64, 0u64), |(sum, n), &x| (sum + x, n + 1))
    };
    let (sum, n) = kept(ruler);
    let mean_kept = sum as f64 / n as f64;
    (0..32)
        .filter_map(|gc| {
            let (sum, n) = kept(&ruler_by_gc[gc]);
            (ruler_by_gc[gc].len() >= 1000 && sum >= 10000).then(|| (gc, sum as f64 / n as f64 / mean_kept))
        })
        .collect()
}
// autosomal denominator classes
fn is_depth_class(name: &str) -> bool {
    name == "chr2depth" || name.starts_with("depth")
}
// sum_gc of each class, and true if the G+C curve was built.
// Array classes: sum of count / curve(G+C), scaled to the G+C mix of the ruler. Ruler and depth classes: raw sum.
fn gc_sums(names: &[String], gc: &[Vec<u8>], counts: &[Vec<u64>]) -> (Vec<f64>, bool) {
    let y = names
        .iter()
        .position(|x| x == "ysinglecopy")
        .unwrap_or_else(|| die("no ysinglecopy class"));
    let raw: Vec<f64> = counts.iter().map(|v| v.iter().sum::<u64>() as f64).collect();
    let mut ruler_by_gc: Vec<Vec<u64>> = vec![vec![]; 32];
    for (&g, &x) in gc[y].iter().zip(&counts[y]) {
        ruler_by_gc[g as usize].push(x);
    }
    let bins = gc_bins(&counts[y], &ruler_by_gc);
    if bins.is_empty() {
        eprintln!("WARNING gc: ysinglecopy sum {}: too few counts for a GC curve, sum_gc = sum; estimate per ruler if y_dose is above the model no-Y cut, else at haploid autosomal depth", raw[y]);
        return (raw, false);
    }
    if bins.iter().any(|&(_, v)| !(v > 0.0)) {
        die("gc: bad curve");
    }
    // a G+C bin without a value takes the nearest bin
    let curve: Vec<f64> = (0..32)
        .map(|g| {
            let mut best = bins[0];
            for &e in &bins[1..] {
                if (e.0 as i64 - g as i64).abs() < (best.0 as i64 - g as i64).abs() {
                    best = e;
                }
            }
            best.1
        })
        .collect();
    // mean curve value over the ruler k-mers
    let mut ruler_mix = 0.0;
    for g in 0..32 {
        ruler_mix += ruler_by_gc[g].len() as f64 * curve[g];
    }
    ruler_mix /= gc[y].len() as f64;
    eprintln!(
        "gc r(GC 10..20) {}",
        curve[10..21]
            .iter()
            .map(|v| format!("{v:.3}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    let sums = names
        .iter()
        .enumerate()
        .map(|(k, name)| {
            if name != "ysinglecopy" && !is_depth_class(name) {
                let mut sum = 0.0;
                for (&g, &x) in gc[k].iter().zip(&counts[k]) {
                    sum += x as f64 / curve[g as usize];
                }
                sum * ruler_mix
            } else {
                raw[k]
            }
        })
        .collect();
    (sums, true)
}

static NAMES: std::sync::OnceLock<Vec<Vec<u8>>> = std::sync::OnceLock::new();
static BINS: AtomicBool = AtomicBool::new(false); // --bins given: aligned input only
static EQB: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0); // bases coded '=' (same as reference): not readable, k-mers across them lost
static READS: std::sync::OnceLock<String> = std::sync::OnceLock::new(); // --reads <file>: write marker reads; off by default, they are read sequences
static QUAL: std::sync::OnceLock<(u8, Vec<i64>)> = std::sync::OnceLock::new(); // --qual: min Q (0 = weight), ln P(base right) x 1e9 per Q
fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
enum Src {
    Stdin,
    Hts {
        input: String,
        reference: Option<String>,
        bed: Option<String>,
    },
}
type Tx = mpsc::SyncSender<(u64, Vec<u8>)>;
// reads into SAM-like chunks cut at line ends
fn spawn_reader(src: Src, tx: Tx, prx: mpsc::Receiver<Vec<u8>>, thr: usize) -> std::thread::JoinHandle<u64> {
    std::thread::Builder::new()
        .name("yr-read".into())
        .spawn(move || match src {
            Src::Stdin => read_stdin(tx, prx),
            Src::Hts { input, reference, bed } => read_hts(&input, reference.as_deref(), bed.as_deref(), thr, tx, prx),
        })
        .unwrap()
}
fn read_stdin(tx: Tx, prx: mpsc::Receiver<Vec<u8>>) -> u64 {
    let mut carry: Vec<u8> = Vec::new();
    let (mut ci, mut bytes) = (0u64, 0u64);
    loop {
        let mut buf = prx
            .try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(chunk() + (1 << 17)));
        buf.clear();
        buf.extend_from_slice(&carry);
        carry.clear();
        let start = buf.len();
        buf.reserve(chunk());
        let bp = unsafe { buf.as_mut_ptr().add(start) }; // read(2) into spare capacity
        let mut got = 0;
        while got < chunk() {
            let k = unsafe { libc::read(0, bp.add(got) as *mut libc::c_void, chunk() - got) };
            if k == 0 {
                break;
            }
            if k < 0 {
                let e = io::Error::last_os_error();
                if e.kind() == io::ErrorKind::Interrupted {
                    continue;
                }
                die(&format!("stdin: {e}"));
            }
            got += k as usize;
        }
        unsafe {
            buf.set_len(start + got);
        }
        bytes += got as u64;
        if got == 0 {
            if !buf.is_empty() {
                tx.send((ci, buf)).unwrap();
            }
            break;
        }
        match buf.iter().rposition(|&b| b == b'\n') {
            Some(p) => {
                carry.extend_from_slice(&buf[p + 1..]);
                buf.truncate(p + 1);
            }
            None => {
                carry = buf;
                continue;
            }
        } // no line end: line longer than the chunk, read on
        tx.send((ci, buf)).unwrap();
        ci += 1;
    }
    bytes
}
// BED regions merged per contig, in header order
fn regions(bed: &str, names: &[Vec<u8>]) -> Vec<(u32, i64, i64)> {
    let mut by: HashMap<Vec<u8>, Vec<(i64, i64)>> = HashMap::new();
    for l in std::fs::read_to_string(bed)
        .unwrap_or_else(|e| die(&format!("{bed}: {e}")))
        .lines()
    {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 3 || f[0].starts_with('#') || f[0] == "*" {
            continue;
        }
        let (a, b) = (
            f[1].parse::<i64>()
                .unwrap_or_else(|_| die(&format!("bed start {}", f[1]))),
            f[2].parse::<i64>()
                .unwrap_or_else(|_| die(&format!("bed end {}", f[2]))),
        );
        by.entry(f[0].as_bytes().to_vec()).or_default().push((a, b));
    }
    let mut out = vec![];
    for (t, n) in names.iter().enumerate() {
        if let Some(v) = by.get_mut(n) {
            v.sort();
            let mut m: Vec<(i64, i64)> = vec![];
            for &(a, b) in v.iter() {
                match m.last_mut() {
                    Some(l) if a <= l.1 => l.1 = l.1.max(b),
                    _ => m.push((a, b)),
                }
            }
            for (a, b) in m {
                out.push((t as u32, a, b));
            }
        }
    }
    let found: usize = by.keys().filter(|k| names.contains(k)).count();
    if found < by.len() {
        eprintln!(
            "WARNING bed: {} of {} contigs not in the input header",
            by.len() - found,
            by.len()
        );
    }
    out
}
// one record: flag u16, tid i32, pos i64, mapq u8, l_seq u32, 4-bit seq; 0xFF opens a binary chunk
fn pack(buf: &mut Vec<u8>, rec: &rust_htslib::bam::Record) {
    if buf.is_empty() {
        buf.push(0xFF);
    }
    let l = rec.seq_len();
    buf.extend_from_slice(&rec.flags().to_le_bytes());
    buf.extend_from_slice(&rec.tid().to_le_bytes());
    buf.extend_from_slice(&rec.pos().to_le_bytes());
    buf.push(rec.mapq());
    buf.extend_from_slice(&(l as u32).to_le_bytes());
    buf.extend_from_slice(&rec.seq().encoded[..(l + 1) / 2]);
    if QUAL.get().is_some() {
        buf.extend_from_slice(rec.qual());
    } // --qual: one phred byte per base
}
// SAM/BAM/CRAM without index, FASTQ/FASTA (plain or gzip), list a,b: one reader per file, start to end
fn read_seq(input: &str, reference: Option<&str>, thr: usize, tx: Tx, prx: mpsc::Receiver<Vec<u8>>) -> u64 {
    use rust_htslib::{bam, htslib as h};
    let cs = |s: &str| std::ffi::CString::new(s).unwrap_or_else(|_| die(&format!("{s}: NUL in the name")));
    let newbuf = || {
        prx.try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(chunk() + (1 << 17)))
    };
    let (mut rec, mut buf, mut ci, mut n) = (bam::Record::new(), newbuf(), 0u64, 0u64);
    buf.clear();
    for f in input.split(',') {
        let fp = unsafe { h::hts_open(cs(f).as_ptr(), c"r".as_ptr()) };
        if fp.is_null() {
            die(&format!("{f}: cannot open"));
        }
        if unsafe { (*fp).format.format } == h::htsExactFormat_empty_format {
            eprintln!("WARNING {f}: empty file, 0 reads");
            unsafe {
                h::hts_close(fp);
            }
            continue;
        }
        let kind = match unsafe { (*fp).format.format } {
            h::htsExactFormat_sam => "SAM",
            h::htsExactFormat_bam => "BAM",
            h::htsExactFormat_cram => "CRAM",
            h::htsExactFormat_fastq_format => "FASTQ",
            h::htsExactFormat_fasta_format => "FASTA",
            _ => die(&format!("{f}: not SAM, BAM, CRAM, FASTQ or FASTA")),
        };
        if BINS.load(Relaxed) && matches!(kind, "FASTQ" | "FASTA") {
            die(&format!(
                "{f}: {kind} has no read positions; --bins needs SAM, BAM or CRAM input"
            ));
        }
        if matches!(kind, "BAM" | "CRAM") {
            eprintln!("WARNING {f}: no index: one reader, start to end (slower)");
        }
        unsafe {
            if let Some(r) = reference.filter(|_| kind == "CRAM") {
                if h::hts_set_fai_filename(fp, cs(r).as_ptr()) != 0 {
                    die(&format!("reference {r}: cannot use it"));
                }
            }
            if h::hts_set_threads(fp, thr.max(1) as i32) != 0 {
                eprintln!("WARNING {f}: no htslib threads");
            }
            if kind == "CRAM"
                && h::hts_set_opt(
                    fp,
                    h::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
                    0x2 | 0x4 | 0x8 | 0x10 | 0x200 | if QUAL.get().is_some() { 0x400 } else { 0 },
                ) != 0
            {
                eprintln!("WARNING {f}: cram required fields not set");
            }
        }
        let hd = unsafe { h::sam_hdr_read(fp) };
        if hd.is_null() {
            die(&format!("{f}: cannot read the header"));
        }
        let names: Vec<Vec<u8>> = (0..unsafe { h::sam_hdr_nref(hd) })
            .map(|t| {
                unsafe { std::ffi::CStr::from_ptr(h::sam_hdr_tid2name(hd, t)) }
                    .to_bytes()
                    .to_vec()
            })
            .collect();
        if NAMES.get_or_init(|| names.clone()) != &names {
            die(&format!("{f}: contigs differ from the first input"));
        }
        let mut m = 0u64;
        loop {
            let r = unsafe { h::sam_read1(fp, hd, rec.inner_mut()) };
            if r == -1 {
                break;
            }
            if r < -1 {
                die(&format!("read {f}: htslib error {r}"));
            }
            if rec.flags() & 0x900 != 0 {
                continue;
            } // secondary/supplementary
            pack(&mut buf, &rec);
            m += 1;
            if buf.len() >= chunk() {
                let mut b = newbuf();
                b.clear();
                tx.send((ci, std::mem::replace(&mut buf, b))).unwrap();
                ci += 1;
            }
        }
        unsafe {
            h::sam_hdr_destroy(hd);
            h::hts_close(fp);
        }
        eprintln!("{f}: {kind}, records read (primary, one reader) {m}");
        n += m;
    }
    if buf.len() > 1 {
        tx.send((ci, buf)).unwrap();
    }
    n
}
// CRAM/BAM: shards (region, previous end on contig) read in parallel, each read once; needs the index
fn read_hts(
    input: &str,
    reference: Option<&str>,
    bed: Option<&str>,
    thr: usize,
    tx: Tx,
    prx: mpsc::Receiver<Vec<u8>>,
) -> u64 {
    use rust_htslib::bam::{self, Read};
    use std::sync::atomic::{AtomicU64, AtomicUsize};
    let idx = || {
        if input.contains("://") {
            bam::IndexedReader::from_url(&input.parse().unwrap_or_else(|e| die(&format!("{input}: {e}"))))
        } else {
            bam::IndexedReader::from_path(input)
        }
    };
    let cfg = |mut r: bam::IndexedReader, t: usize| {
        if let Some(f) = reference {
            r.set_reference(f)
                .unwrap_or_else(|e| die(&format!("reference {f}: {e}")));
        }
        r.set_threads(t.max(1))
            .unwrap_or_else(|e| die(&format!("threads: {e}")));
        if let Err(e) = r.set_cram_options(
            rust_htslib::htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            0x2 | 0x4 | 0x8 | 0x10 | 0x200 | if QUAL.get().is_some() { 0x400 } else { 0 },
        ) {
            eprintln!("WARNING cram required fields: {e}");
        }
        r
    };
    let open = |t: usize| cfg(idx().unwrap_or_else(|e| die(&format!("{input}: {e}"))), t);
    let h = match idx() {
        Ok(r) => cfg(r, 1),
        Err(_) if bed.is_none() => return read_seq(input, reference, thr, tx, prx),
        Err(e) => die(&format!("{input}: {e} (-b needs an indexed BAM or CRAM)")),
    };
    let names: Vec<Vec<u8>> = h.header().target_names().iter().map(|n| n.to_vec()).collect();
    let _ = NAMES.set(names.clone());
    let rg: Vec<(u32, i64, i64)> = match bed {
        Some(b) => regions(b, &names),
        None => (0..names.len() as u32)
            .flat_map(|t| {
                let l = h.header().target_len(t).unwrap_or(0) as i64;
                (0..l.max(1))
                    .step_by(10_000_000)
                    .map(move |a| (t, a, (a + 10_000_000).min(l)))
            })
            .collect(),
    };
    let mut sh: Vec<Option<(u32, i64, i64, i64)>> = rg
        .iter()
        .enumerate()
        .map(|(i, &(t, a, e))| {
            Some((
                t,
                a,
                e,
                if i > 0 && rg[i - 1].0 == t {
                    rg[i - 1].2
                } else {
                    i64::MIN
                },
            ))
        })
        .collect();
    sh.push(None); // unmapped
    let nr = (thr / 2).max(1).min(sh.len());
    eprintln!(
        "shards {} ({}) + unmapped, readers {nr}, htslib threads each {}",
        sh.len() - 1,
        if bed.is_some() { "bed regions" } else { "10 Mb windows" },
        (thr / nr).max(1)
    );
    let (next, ci, nrec, prx) = (
        AtomicUsize::new(0),
        AtomicU64::new(0),
        AtomicU64::new(0),
        Mutex::new(prx),
    );
    std::thread::scope(|sc| {
        for _ in 0..nr {
            sc.spawn(|| {
                let mut r = open(thr / nr);
                let mut rec = bam::Record::new();
                let mut n = 0u64;
                let newbuf = || {
                    prx.lock()
                        .unwrap()
                        .try_recv()
                        .unwrap_or_else(|_| Vec::with_capacity(chunk() + (1 << 17)))
                };
                let mut buf = newbuf();
                buf.clear();
                loop {
                    let i = next.fetch_add(1, Relaxed);
                    if i >= sh.len() {
                        break;
                    }
                    let (t, a, pe) = match sh[i] {
                        Some((t, a, _, pe)) => (t as i32, a, pe),
                        None => (-1, i64::MIN, i64::MIN),
                    };
                    for k in 1..=6 {
                        // a failed fetch reads no record: reopen and retry (remote files)
                        let x = match sh[i] {
                            Some((t, a, e, _)) => r.fetch((t as i32, a, e)),
                            None => r.fetch(bam::FetchDefinition::Unmapped),
                        };
                        let Err(x) = x else { break };
                        let w = match sh[i] {
                            Some((t, a, e, _)) => format!("{}:{a}-{e}", String::from_utf8_lossy(&names[t as usize])),
                            None => "unmapped".into(),
                        };
                        if k == 6 {
                            die(&format!("fetch {w}: {x}"));
                        }
                        eprintln!("WARNING fetch {w}: {x}; retry {k} of 5 in {} s", 10 * k);
                        std::thread::sleep(std::time::Duration::from_secs(10 * k));
                        r = open(thr / nr);
                    }
                    while let Some(x) = r.read(&mut rec) {
                        x.unwrap_or_else(|e| die(&format!("read {input}: {e}")));
                        let (fl, p) = (rec.flags(), rec.pos());
                        if fl & 0x900 != 0 || (t >= 0 && p < a && p < pe) {
                            continue;
                        } // secondary/supplementary, or seen in the previous shard
                        pack(&mut buf, &rec);
                        n += 1;
                        if buf.len() >= chunk() {
                            let mut b = newbuf();
                            b.clear();
                            tx.send((ci.fetch_add(1, Relaxed), std::mem::replace(&mut buf, b)))
                                .unwrap();
                        }
                    }
                }
                if buf.len() > 1 {
                    tx.send((ci.fetch_add(1, Relaxed), buf)).unwrap();
                }
                nrec.fetch_add(n, Relaxed);
            });
        }
    });
    let n = nrec.load(Relaxed);
    eprintln!("records read (primary) {n}");
    n
}

fn count(src: Src, fa_path: &str, name: &str, bin: u64, wbins: bool, threads: usize, out: &str) {
    let t0 = Instant::now();
    // reader: stdin read during the marker load (chunks queued, size printed at start), chunks cut at line ends
    let queue: usize = 32;
    eprintln!("queue {queue} chunks of {} MB", chunk() >> 20);
    let (tx, rx) = mpsc::sync_channel::<(u64, Vec<u8>)>(queue);
    let rx = Arc::new(Mutex::new(rx));
    let (ptx, prx) = mpsc::channel::<Vec<u8>>(); // used chunk buffers back to the reader (no new pages per chunk)
    let reader = spawn_reader(src, tx, prx, threads);
    let mut markers = load(fa_path, threads);
    let ncl = markers.names.len();
    let n = markers.key.len();
    let per: Vec<usize> = (0..ncl)
        .map(|j| markers.ocls.iter().filter(|&&c| c as usize == j).count())
        .collect();
    let tab = std::mem::replace(
        &mut markers.tab,
        Table {
            s: Vec::new(),
            cls: Vec::new(),
            mask: 0,
        },
    );
    eprintln!(
        "markers lines {} distinct {} classes {:?} per class {:?}; table {} MB; load {:.1} s; {}",
        markers.order.len(),
        n,
        markers.names,
        per,
        (tab.s.len() * 16) >> 20,
        t0.elapsed().as_secs_f64(),
        rss()
    );
    let bpk: usize = 16; // Bloom bits per key
    let bloom = {
        let mut b = Bloom::new(n, bpk.max(1));
        for &x in &markers.key {
            b.insert(kmer_hash(x));
        }
        Arc::new(b)
    };
    eprintln!("bloom {} MB ({} bits/key)", bloom.b.len() * 8 >> 20, bpk.max(1));
    let tab = Arc::new(tab);
    let cb: u32 = 15; // per-thread cache: 2^cb entries
    let (ftx, frx) = mpsc::sync_channel::<(u64, Vec<u8>)>(threads * 4);
    std::fs::create_dir_all(out).unwrap_or_else(|e| die(&format!("{out}: {e}")));
    let any_path = READS.get().cloned();
    let writer = std::thread::Builder::new()
        .name("yr-write".into())
        .spawn(move || {
            let Some(p) = any_path else {
                for _ in frx {}
                return;
            }; // no read sequences on disk by default
            let mut w = BufWriter::with_capacity(
                1 << 22,
                std::fs::File::create(&p).unwrap_or_else(|e| die(&format!("{p}: {e}"))),
            );
            let mut next = 0u64;
            let mut held: std::collections::BTreeMap<u64, Vec<u8>> = std::collections::BTreeMap::new(); // marker reads in input order
            for (i, b) in frx {
                held.insert(i, b);
                while let Some(b) = held.remove(&next) {
                    w.write_all(&b).unwrap_or_else(|e| die(&format!("write {p}: {e}")));
                    next += 1;
                }
            }
            if !held.is_empty() {
                die("marker reads: chunks missing");
            }
            w.flush().unwrap();
        })
        .unwrap();
    let workers: Vec<_> = (0..threads)
        .map(|_| {
            let (rx, ftx, ptx, tab, bloom) = (rx.clone(), ftx.clone(), ptx.clone(), tab.clone(), bloom.clone());
            std::thread::Builder::new().name("yr-count".into()).spawn(move || {
                let mut loc = Local {
                    contigs: HashMap::new(),
                    cname: vec![],
                    bins: HashMap::new(),
                    rk: vec![],
                    rv: vec![],
                    reads: 0,
                    kmers: 0,
                    hits: 0,
                    cmiss: 0,
                    bpos: 0,
                    cache: vec![
                        Ent {
                            key: EMPTY,
                            sc: 0,
                            pend: 0,
                            cls: 0
                        };
                        1 << cb
                    ],
                };
                let mut kb: Vec<(u64, u64, bool)> = Vec::with_capacity(1024);
                let mut sq: Vec<u8> = Vec::with_capacity(1 << 16);
                let mut cx = Ctx {
                    bloom: &bloom,
                    tab: &tab,
                    ncl,
                    bin,
                    h: vec![0; ncl],
                    last: (Vec::new(), u64::MAX, 0),
                };
                loop {
                    let c = {
                        let r = rx.lock().unwrap();
                        r.recv()
                    };
                    let Ok((ci, c)) = c else { break };
                    let mut fa = Vec::with_capacity(1 << 16);
                    if c.first() == Some(&0xFF) {
                        work_bin(&c, &mut loc, &mut cx, &mut fa, &mut kb, &mut sq);
                    } else {
                        work(&c, &mut loc, &mut cx, &mut fa, &mut kb);
                    }
                    ftx.send((ci, fa)).unwrap();
                    let _ = ptx.send(c); // fails only after the reader ended: buffer then freed
                }
                for e in &loc.cache {
                    evict(e, &tab);
                }
                loc.bins = HashMap::new();
                loc.cache = vec![];
                loc // lookup maps not needed after counting
            })
        })
        .collect();
    drop(ftx);
    let bytes_in = reader.join().unwrap();
    let locs: Vec<Local> = workers.into_iter().map(|w| w.unwrap().join().unwrap()).collect();
    writer.join().unwrap();
    let t_count = t0.elapsed().as_secs_f64();
    eprintln!("counting done; {}", rss());
    let tp = Instant::now();
    let lap = |w: &str| eprintln!("post {w} at {:.2} s", tp.elapsed().as_secs_f64());
    // bins: global contig ids, rows sorted by (contig name, bin); worker rows summed while writing (no merged copy)
    let w = 3 + 3 * ncl;
    let (mut contig_names, mut contig_ids): (Vec<Vec<u8>>, HashMap<Vec<u8>, u32>) = (vec![], HashMap::new());
    let (mut reads, mut kmers, mut hits, mut cmiss, mut bpos) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut bin_rows: Vec<(u32, u64, u32, u32)> = Vec::with_capacity(locs.iter().map(|l| l.rk.len()).sum()); // (contig, bin, worker, row)
    for (wi, l) in locs.iter().enumerate() {
        reads += l.reads;
        kmers += l.kmers;
        cmiss += l.cmiss;
        bpos += l.bpos;
        hits += l.hits;
        let g: Vec<u32> = l
            .cname
            .iter()
            .map(|c| {
                *contig_ids.entry(c.clone()).or_insert_with(|| {
                    contig_names.push(c.clone());
                    contig_names.len() as u32 - 1
                })
            })
            .collect();
        for (r, &(ci, b)) in l.rk.iter().enumerate() {
            bin_rows.push((g[ci as usize], b, wi as u32, r as u32));
        }
    }
    if reads == 0 {
        die("0 reads");
    }
    let mut rank = vec![0u32; contig_names.len()];
    {
        let mut o: Vec<usize> = (0..contig_names.len()).collect();
        o.sort_by(|&a, &b| contig_names[a].cmp(&contig_names[b]));
        for (r, &i) in o.iter().enumerate() {
            rank[i] = r as u32;
        }
    }
    bin_rows.sort_unstable_by_key(|&(c, b, wi, r)| (rank[c as usize], b, wi, r));
    let n_bins = bin_rows.len()
        - bin_rows
            .windows(2)
            .filter(|p| p[0].0 == p[1].0 && p[0].1 == p[1].1)
            .count();
    eprintln!(
        "bin records over workers {} (merged {n_bins}); {}",
        bin_rows.len(),
        rss()
    );
    eprintln!("reads {reads} bins {n_bins} bytes {} k-mers {kmers} cache misses {cmiss} ({:.2}%) bloom positives {bpos} hits {hits}; counting {:.1} s", bytes_in, 100.0 * cmiss as f64 / kmers.max(1) as f64, t_count);
    lap("bins merged");
    let eq = EQB.load(Relaxed);
    if eq > 0 {
        die(&format!("{eq} read bases are written as '=' (same as the reference); yruler reads bases, not alignments, so it cannot recover them and the counts would be low. Write the file with explicit bases and run again"));
    } // no output with unreadable bases
    let mut kmer_counts: Vec<u64> = vec![0; n];
    for sl in tab.s.iter() {
        if sl.key != EMPTY {
            kmer_counts[sl.idx as usize] = sl.cnt.load(Relaxed) as u64;
        }
    }
    drop(tab);
    drop(bloom); // counts extracted: table and filter freed
                 // per class lists in marker order
    let mut counts_by_class: Vec<Vec<u64>> = vec![vec![]; ncl];
    let mut gc_by_class: Vec<Vec<u8>> = vec![vec![]; ncl];
    for (k, &i) in markers.order.iter().enumerate() {
        let c = markers.ocls[k] as usize;
        counts_by_class[c].push(kmer_counts[i as usize]);
        gc_by_class[c].push(markers.gc[k]);
    }
    // check: bins hits = counts per class
    for j in 0..ncl {
        let s: u64 = counts_by_class[j].iter().sum();
        let bh: u64 = locs
            .iter()
            .map(|l| l.rv.iter().skip(3 + 3 * j).step_by(w).sum::<u64>())
            .sum();
        eprintln!("check {} counts {s} bins {bh}", markers.names[j]);
        if s != bh {
            die(&format!("{} bins hits {bh} vs counts {s}", markers.names[j]));
        }
    }
    let wr = |p: &str| {
        BufWriter::with_capacity(
            1 << 22,
            std::fs::File::create(format!("{out}/{p}")).unwrap_or_else(|e| die(&format!("{p}: {e}"))),
        )
    };
    let t1 = Instant::now();
    {
        // counts.tsv: rounds of thr parallel blocks (1M lines each), written in order
        let (key, kmer_counts, per) = (&markers.key, &kmer_counts, 1usize << 20);
        let mut o = wr("counts.tsv");
        for r0 in (0..n).step_by(per * threads) {
            let parts: Vec<Vec<u8>> = std::thread::scope(|sc| {
                let hs: Vec<_> = (r0..(r0 + per * threads).min(n))
                    .step_by(per)
                    .map(|a| {
                        sc.spawn(move || {
                            let mut b = Vec::with_capacity(per * 40);
                            for i in a..(a + per).min(n) {
                                b.extend_from_slice(&decode(key[i]));
                                b.push(b'\t');
                                itoa(&mut b, kmer_counts[i]);
                                b.push(b'\n');
                            }
                            b
                        })
                    })
                    .collect();
                hs.into_iter().map(|h| h.join().unwrap()).collect()
            });
            for b in &parts {
                o.write_all(b)
                    .unwrap_or_else(|e| die(&format!("write counts.tsv: {e}")));
            }
        }
        o.flush().unwrap();
    }
    eprintln!("write counts {:.2} s", t1.elapsed().as_secs_f64());
    if wbins {
        let mut o = wr("bins.tsv");
        let mut v = vec![0u64; w];
        let mut b: Vec<u8> = Vec::with_capacity(1 << 20);
        write!(o, "contig\tstart\treads\treads_mq20\tduplicates").unwrap();
        for c in &markers.names {
            write!(o, "\t{c}_hits\t{c}_hits_mq20\t{c}_reads").unwrap();
        }
        writeln!(o).unwrap();
        for (i, &(c, bn, wi, r)) in bin_rows.iter().enumerate() {
            let x = &locs[wi as usize].rv[r as usize * w..(r as usize + 1) * w];
            for (a, y) in v.iter_mut().zip(x) {
                *a += y;
            }
            if bin_rows.get(i + 1).map_or(true, |n| n.0 != c || n.1 != bn) {
                b.extend_from_slice(&contig_names[c as usize]);
                b.push(b'\t');
                itoa(&mut b, bn * bin);
                for a in v.iter_mut() {
                    b.push(b'\t');
                    itoa(&mut b, *a);
                    *a = 0;
                }
                b.push(b'\n');
                if b.len() >= 1 << 20 {
                    o.write_all(&b).unwrap_or_else(|e| die(&format!("write bins.tsv: {e}")));
                    b.clear();
                }
            }
        }
        o.write_all(&b).unwrap_or_else(|e| die(&format!("write bins.tsv: {e}")));
        o.flush().unwrap();
    }
    if !markers.names.iter().any(|x| x == "ysinglecopy") {
        die("missing class ysinglecopy");
    }
    let mut depth_means: Vec<f64> = markers
        .names
        .iter()
        .enumerate()
        .filter(|(_, n)| is_depth_class(n))
        .map(|(j, _)| counts_by_class[j].iter().sum::<u64>() as f64 / counts_by_class[j].len() as f64)
        .collect();
    if depth_means.is_empty() {
        die("missing denominator class (chr2depth or depth*)");
    }
    depth_means.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ix = |c: &str| markers.names.iter().position(|x| x == c).unwrap();
    let ruler_mean =
        counts_by_class[ix("ysinglecopy")].iter().sum::<u64>() as f64 / counts_by_class[ix("ysinglecopy")].len() as f64;
    let depth_mean = if depth_means.len() % 2 == 1 {
        depth_means[depth_means.len() / 2]
    } else {
        (depth_means[depth_means.len() / 2 - 1] + depth_means[depth_means.len() / 2]) / 2.0
    }; // median of the window means
    if depth_mean == 0.0 {
        die("denominator mean 0");
    }
    lap("counts and bins written");
    let (sum_gc, gc_curve) = gc_sums(&markers.names, &gc_by_class, &counts_by_class);
    lap("gc");
    {
        let mut o = wr("summary.tsv");
        writeln!(
            o,
            "#y_dose={}\tysinglecopy_depth={}\tdepth_mean={}\tgc_curve={}",
            py_float(2.0 * ruler_mean / depth_mean),
            py_float(ruler_mean),
            py_float(depth_mean),
            gc_curve as u8
        )
        .unwrap();
        writeln!(o, "sample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc").unwrap();
        for (j, c) in markers.names.iter().enumerate() {
            let v = &counts_by_class[j];
            let (_, ms) = median(v);
            let z = v.iter().filter(|&&x| x == 0).count();
            writeln!(
                o,
                "{name}\t{reads}\t{c}\t{}\t{}\t{ms}\t{:.4}\t{}",
                v.len(),
                v.iter().sum::<u64>(),
                z as f64 / v.len() as f64,
                py_float(sum_gc[j])
            )
            .unwrap();
        } // full sum_gc: estimate reads it back
        o.flush().unwrap();
    }
    eprintln!(
        "y_dose {:?} ysinglecopy_depth {:?} depth_mean {:?} (median of {} denominator classes); total {:.1} s; {}",
        2.0 * ruler_mean / depth_mean,
        ruler_mean,
        depth_mean,
        depth_means.len(),
        t0.elapsed().as_secs_f64(),
        rss()
    );
}

const TECH: [&str; 4] = ["illumina", "hifi", "ont", "hic"];
// estimate_bp = K x relative, with K of the class for technology x (model column 3, empty = illumina).
// y_found (scaled y_dose above the model no_y_dose): relative = sum_gc(class) / sum(ysinglecopy).
// Without a G+C curve, sum_gc = sum. no_y (y_dose at or below no_y_dose, or no no_y_dose line and gc_curve=0):
// relative = sum(class) / (n_kmers(ysinglecopy) x depth_mean / 2), the signal at haploid autosomal depth.
// depth_mean: median over the denominator classes (chr2depth or depth*) of their mean count.
// y_dose = 2 x ysinglecopy mean / depth_mean, divided in estimates.tsv by the model y_dose_scale.
// Model column 4 (optional): the reported value as a signed sum of classes (e.g. dyz1+dyz2+dyz18s+yq12rest2), none = not reported.
fn estimate(summary: &str, model: &str, x: &str, out: &mut dyn Write) {
    let read = |p: &str| std::fs::read_to_string(p).unwrap_or_else(|e| die(&format!("{p}: {e}")));
    let (summary_text, model_text) = (read(summary), read(model));

    // summary.tsv: header line, then one row per class
    let mut y_dose_text = String::new();
    let mut gc_curve = true;
    let mut depth_mean = f64::NAN;
    let mut ruler_sum = f64::NAN;
    let mut ruler_kmers = f64::NAN;
    let mut rows: Vec<(String, String, f64)> = vec![]; // (sample, class, sum_gc)
    let mut raw_sum: HashMap<String, f64> = HashMap::new();
    for line in summary_text.lines() {
        if let Some(header) = line.strip_prefix('#') {
            let get = |key: &str| header.split('\t').find_map(|f| f.strip_prefix(key)).unwrap_or("");
            y_dose_text = get("y_dose=").to_string();
            gc_curve = get("gc_curve=") != "0";
            depth_mean = [get("depth_mean="), get("chr2depth_mean=")]
                .iter()
                .find(|v| !v.is_empty())
                .map_or(f64::NAN, |v| v.parse().unwrap_or(f64::NAN));
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 || f[0] == "sample" {
            continue;
        }
        let (sample, class, n_kmers, sum, sum_gc) = (f[0], f[2], f[3], f[4], f[7]);
        let sum_gc: f64 = sum_gc.parse().unwrap_or_else(|_| die(&format!("sum_gc {sum_gc}")));
        if class == "ysinglecopy" {
            ruler_sum = sum.parse().unwrap_or_else(|_| die("ysinglecopy sum"));
            ruler_kmers = n_kmers.parse().unwrap_or(f64::NAN);
        }
        raw_sum.insert(class.to_string(), sum.parse::<f64>().unwrap_or(f64::NAN));
        rows.push((sample.to_string(), class.to_string(), sum_gc));
    }

    // model.tsv: header values
    let model_value = |key: &str| {
        model_text.lines().filter_map(|l| l.strip_prefix(key)).next().map(|v| {
            v.trim()
                .parse::<f64>()
                .unwrap_or_else(|_| die(&format!("model {} {v}", &key[1..key.len() - 1])))
        })
    };
    let y_dose_scale = model_value("#y_dose_scale=").unwrap_or(1.0); // y_dose of the reference male in the model cohort
    let no_y_dose = model_value("#no_y_dose="); // highest female y_dose in the model cohort

    // Y found: y_dose above the no-Y cut, both with 4 decimals
    let y_found = match (no_y_dose, y_dose_text.parse::<f64>()) {
        (Some(cut), Ok(d)) => (d / y_dose_scale * 1e4).round() / 1e4 > cut,
        _ => gc_curve,
    };
    if y_found && !(ruler_sum > 0.0) {
        die("no ysinglecopy sum > 0 in summary");
    }
    if !y_found && !(ruler_kmers > 0.0 && depth_mean > 0.0) {
        die("the no_y estimate needs ysinglecopy n_kmers and depth_mean in the summary");
    }
    if y_found && !gc_curve {
        eprintln!("WARNING no G+C curve: Y found (y_dose above the no-Y cut), estimates without G+C correction");
    }
    let y_dose = y_dose_text
        .parse::<f64>()
        .map_or(y_dose_text.clone(), |d| format!("{:.4}", d / y_dose_scale));
    let y_status = if y_found { "y_found" } else { "no_y" };

    // model.tsv: rows of technology x
    fn tech_of<'a>(f: &[&'a str]) -> &'a str {
        f.get(2).copied().filter(|t| !t.is_empty()).unwrap_or("illumina")
    }
    let model_rows: Vec<Vec<&str>> = model_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(|t| t.trim()).collect::<Vec<&str>>())
        .filter(|f| f.len() >= 2 && f[0] != "class")
        .collect();
    let mut k_of: HashMap<String, f64> = HashMap::new();
    let mut reported: Vec<(String, String)> = vec![]; // (class, signed sum of classes)
    for f in model_rows.iter().filter(|f| tech_of(f) == x) {
        k_of.insert(
            f[0].to_string(),
            f[1].parse().unwrap_or_else(|_| die(&format!("model K {}", f[1]))),
        );
        let terms = f.get(3).copied().filter(|t| !t.is_empty()).unwrap_or(f[0]);
        if terms != "none" {
            reported.push((f[0].to_string(), terms.replace('-', "+-")));
        }
    }
    if k_of.is_empty() {
        let mut techs: Vec<&str> = model_rows.iter().map(|f| tech_of(f)).collect();
        techs.sort();
        techs.dedup();
        die(&format!(
            "model {model} has no K for technology {x}: it is calibrated for {}; use a model calibrated for {x}",
            techs.join(", ")
        ));
    }
    for c in k_of.keys() {
        if !rows.iter().any(|r| &r.1 == c) {
            die(&format!("model class {c} missing in summary"));
        }
    }
    if x == "hic" {
        eprintln!("WARNING technology {x} is experimental");
    }

    // one output row per reported class
    let by_class: HashMap<&str, (&str, f64)> = rows.iter().map(|(s, c, v)| (c.as_str(), (s.as_str(), *v))).collect();
    let k = |c: &str| {
        *k_of
            .get(c)
            .unwrap_or_else(|| die(&format!("length term {c}: no K for technology {x}")))
    };
    let relative = |c: &str| {
        if y_found {
            by_class[c].1 / ruler_sum
        } else {
            raw_sum[c] / (ruler_kmers * depth_mean / 2.0)
        }
    };
    let estimate_bp = |c: &str| k(c) * relative(c);
    writeln!(
        out,
        "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status"
    )
    .unwrap();
    for (class, terms) in &reported {
        let bp: f64 = terms
            .split('+')
            .filter(|t| !t.is_empty())
            .map(|t| {
                t.strip_prefix('-')
                    .map_or_else(|| estimate_bp(t), |neg| -estimate_bp(neg))
            })
            .sum();
        let (sample, sum_gc) = by_class[class.as_str()];
        writeln!(
            out,
            "{sample}\t{class}\t{bp:.0}\t{}\t{sum_gc:.1}\t{ruler_sum:.0}\t{y_dose}\t{y_status}",
            py_float(relative(class))
        )
        .unwrap();
    }
    eprintln!(
        "estimate technology {x} classes {} y_dose {y_dose} (scale {y_dose_scale}) y_status {y_status}",
        reported.len()
    );
}

const USAGE: &str = "usage:
  yruler -i <cram|bam|sam|fastq|fasta|-> -m <markers> -o <outdir> [-r <ref.fa>] [-b <regions.bed>] [-M <model.tsv>] [-x <tech>] [-n <name>] [-t <threads>] [--bins <bp>] [--reads <file>]
  yruler estimate <summary.tsv> <model.tsv> [-x <tech>]
-i -: SAM on stdin; -i a,b: files read one after the other (FASTQ: plain or gzip, interleaved pairs in one file or R1,R2)
-x: illumina (default), hifi, ont, hic (experimental); the model needs K values for that technology (the released model.tsv: illumina)
--bins <bp>: also write bins.tsv, reads and marker hits per <bp> of the reference (aligned input)
--reads <file>: also write the sequences of the reads with a marker hit (off by default)";
fn tech(t: &str) -> &'static str {
    TECH.iter()
        .find(|&&x| x == t)
        .copied()
        .unwrap_or_else(|| die(&format!("-x {t}: use one of {}", TECH.join(", "))))
}
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    // static libcurl: point it to the system CA bundle
    if std::env::var_os("CURL_CA_BUNDLE").is_none() {
        if let Some(p) = [
            "/etc/ssl/certs/ca-certificates.crt",
            "/etc/pki/tls/certs/ca-bundle.crt",
            "/etc/ssl/cert.pem",
        ]
        .iter()
        .find(|p| std::path::Path::new(p).exists())
        {
            unsafe { std::env::set_var("CURL_CA_BUNDLE", p) };
        }
    }
    match args.first().map(|s| s.as_str()) {
        Some("estimate") => {
            let technology = match args.len() {
                3 => "illumina",
                5 if args[3] == "-x" => tech(&args[4]),
                _ => die(USAGE),
            };
            estimate(&args[1], &args[2], technology, &mut io::stdout().lock());
            return;
        }
        Some("--version") | Some("-V") => {
            println!("yruler {}", env!("CARGO_PKG_VERSION"));
            return;
        }
        None | Some("-h") | Some("--help") => {
            eprintln!("{USAGE}");
            std::process::exit(if args.is_empty() { 1 } else { 0 });
        }
        _ => {}
    }
    let mut opts: HashMap<&str, String> = HashMap::new();
    let mut i = 0;
    while i < args.len() {
        let k = args[i].as_str();
        if !matches!(
            k,
            "-i" | "-m" | "-o" | "-r" | "-b" | "-M" | "-x" | "-n" | "-t" | "--bins" | "--reads"
        ) {
            die(&format!("unknown option {k}\n{USAGE}"));
        }
        let v = args.get(i + 1).unwrap_or_else(|| die(&format!("{k} needs a value")));
        if opts.insert(k, v.clone()).is_some() {
            die(&format!("{k} given twice (several inputs: -i a,b)"));
        }
        i += 2;
    }
    let need = |k: &str| {
        opts.get(k)
            .cloned()
            .unwrap_or_else(|| die(&format!("missing {k}\n{USAGE}")))
    };
    let (input, markers, out) = (need("-i"), need("-m"), need("-o"));
    let technology = opts.get("-x").map_or("illumina", |t| tech(t));
    let threads = opts
        .get("-t")
        .map_or(std::thread::available_parallelism().map_or(1, |n| n.get()), |t| {
            t.parse().unwrap_or_else(|_| die("-t"))
        });
    if let Some(p) = opts.get("--reads") {
        READS.set(p.clone()).unwrap();
    }
    let bins: Option<u64> = opts.get("--bins").map(|b| {
        b.parse()
            .ok()
            .filter(|&n| n > 0)
            .unwrap_or_else(|| die("--bins: bin size in bp"))
    });
    let name = opts.get("-n").cloned().unwrap_or_else(|| {
        let b = std::path::Path::new(&input)
            .file_name()
            .map_or("sample".into(), |f| f.to_string_lossy().to_string());
        b.split('.').next().unwrap_or("sample").to_string()
    });
    let src = if input == "-" {
        Src::Stdin
    } else {
        Src::Hts {
            input: input.clone(),
            reference: opts.get("-r").cloned(),
            bed: opts.get("-b").cloned(),
        }
    };
    eprintln!(
        "yruler {} input {input} markers {markers} name {name} threads {threads} technology {technology} regions {}",
        env!("CARGO_PKG_VERSION"),
        opts.get("-b").map_or("whole file", |s| s.as_str())
    );
    // --qual hidden: not calibrated with the released model; future work: test if base qualities help
    // if let Some(q) = o.get("--qual") { let m: u8 = if q == "weight" { 0 } else { q.strip_prefix("min:").and_then(|v| v.parse().ok()).filter(|&v: &u8| v > 0).unwrap_or_else(|| die("--qual: min:<Q> or weight")) };
    // let _ = QUAL.set((m, (0..256).map(|q| if q == 255 { 0 } else if q == 0 { -1_000_000_000_000_000 } else { ((1.0 - 10f64.powf(-(q as f64) / 10.0)).ln() * 1e9) as i64 }).collect())); eprintln!("base qualities: {}", if m == 0 { "k-mers kept with P(all bases right)".to_string() } else { format!("k-mers with a base below Q{m} dropped") }); }
    BINS.store(bins.is_some(), Relaxed);
    count(src, &markers, &name, bins.unwrap_or(500), bins.is_some(), threads, &out);
    if let Some(m) = opts.get("-M") {
        let p = format!("{out}/estimates.tsv");
        let mut w = BufWriter::new(std::fs::File::create(&p).unwrap_or_else(|e| die(&format!("{p}: {e}"))));
        estimate(&format!("{out}/summary.tsv"), m, technology, &mut w);
        w.flush().unwrap();
        print!("{}", std::fs::read_to_string(&p).unwrap());
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn masked_superstring() {
        // uppercase starts a k-mer; header name before '_'
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGTa";
        let mut ms = b">c1_x\n".to_vec();
        ms.extend(s.iter().enumerate().map(|(i, b)| {
            if i < 2 {
                b.to_ascii_uppercase()
            } else {
                b.to_ascii_lowercase()
            }
        }));
        ms.push(b'\n');
        let r = super::ms_records(&ms);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0, b"c1");
        let st: Vec<usize> = (0..r[0].1.len() + 1 - super::K)
            .filter(|&i| r[0].1[i].is_ascii_uppercase())
            .collect();
        assert_eq!(st, vec![0, 1]);
        assert_eq!(
            super::canon(&r[0].1[1..32]),
            super::canon(b"CGTACGTACGTACGTACGTACGTACGTACGT")
        );
    }
    #[test]
    fn median_select() {
        // O(n) median = sort median
        let mut x = 7u64;
        for n in 1..400 {
            let v: Vec<u64> = (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (x >> 59) % 13
                })
                .collect();
            let mut s = v.clone();
            s.sort();
            let m = if n % 2 == 1 {
                s[n / 2] as f64
            } else {
                (s[n / 2 - 1] + s[n / 2]) as f64 / 2.0
            };
            assert_eq!(super::median(&v).0, m);
        }
    }
    #[test]
    fn estimate_terms() {
        // signed sum of class lengths; 'none' not reported
        let d = std::env::temp_dir().join(format!("yruler_est_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let (s, m) = (d.join("s.tsv"), d.join("m.tsv"));
        std::fs::write(&s, "#y_dose=1.0\nsample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc\nx\t1\ta\t1\t1\t1\t0\t10.0\nx\t1\tb\t1\t1\t1\t0\t4.0\nx\t1\tysinglecopy\t1\t100\t1\t0\t100.0\n").unwrap();
        std::fs::write(
            &m,
            "class\tK\ttech\tlength\na\t100\tillumina\ta-b\nb\t50\tillumina\tnone\nc\t1\thifi\n",
        )
        .unwrap();
        let mut o = vec![];
        super::estimate(s.to_str().unwrap(), m.to_str().unwrap(), "illumina", &mut o);
        std::fs::remove_dir_all(&d).unwrap();
        assert_eq!(String::from_utf8(o).unwrap(), "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status\nx\ta\t8\t0.1\t10.0\t100\t1.0000\ty_found\n");
    }
    #[test]
    fn estimate_haploid() {
        // no G+C curve: K * sum / (ruler k-mers * chr2depth mean / 2)
        let d = std::env::temp_dir().join(format!("yruler_hap_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let (s, m) = (d.join("s.tsv"), d.join("m.tsv"));
        std::fs::write(&s, "#y_dose=0.0\tysinglecopy_depth=0.0\tchr2depth_mean=20\tgc_curve=0\nsample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc\nx\t1\ta\t1\t3\t1\t0\t3.0\nx\t1\tysinglecopy\t10\t0\t0\t1\t0.0\n").unwrap();
        std::fs::write(&m, "class\tK\ttech\tlength\na\t100\tillumina\ta\n").unwrap();
        let mut o = vec![];
        super::estimate(s.to_str().unwrap(), m.to_str().unwrap(), "illumina", &mut o);
        std::fs::remove_dir_all(&d).unwrap();
        assert_eq!(String::from_utf8(o).unwrap(), "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status\nx\ta\t3\t0.03\t3.0\t0\t0.0000\tno_y\n");
        // 100 * 3 / (10 * 20 / 2)
    }
    #[test]
    fn estimate_no_curve_ruler() {
        // no G+C curve, y_dose above the no-Y cut: ruler basis, sum_gc = sum; below: haploid
        let d = std::env::temp_dir().join(format!("yruler_ncr_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let (s, m) = (d.join("s.tsv"), d.join("m.tsv"));
        let run = |yd: &str| {
            std::fs::write(&s, format!("#y_dose={yd}\tysinglecopy_depth=0.1\tdepth_mean=20\tgc_curve=0\nsample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc\nx\t1\ta\t1\t3\t1\t0\t3.0\nx\t1\tysinglecopy\t10\t30\t0\t1\t30.0\n")).unwrap();
            std::fs::write(
                &m,
                "#y_dose_scale=0.5\n#no_y_dose=0.01\nclass\tK\ttech\tlength\na\t100\tillumina\ta\n",
            )
            .unwrap();
            let mut o = vec![];
            super::estimate(s.to_str().unwrap(), m.to_str().unwrap(), "illumina", &mut o);
            String::from_utf8(o).unwrap()
        };
        let (a, b) = (run("0.5"), run("0.004"));
        std::fs::remove_dir_all(&d).unwrap();
        assert_eq!(a, "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status\nx\ta\t10\t0.1\t3.0\t30\t1.0000\ty_found\n"); // 100 * 3 / 30
        assert_eq!(b, "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status\nx\ta\t3\t0.03\t3.0\t30\t0.0080\tno_y\n");
        // 100 * 3 / (10 * 20 / 2)
    }
    #[test]
    fn estimate_scale_tech() {
        // empty tech = illumina; y_dose / model y_dose_scale
        let d = std::env::temp_dir().join(format!("yruler_sc_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let (s, m) = (d.join("s.tsv"), d.join("m.tsv"));
        std::fs::write(&s, "#y_dose=1.0\nsample\tnreads\tclass\tn_kmers\tsum\tmedian\tfrac_zero\tsum_gc\nx\t1\ta\t1\t1\t1\t0\t10.0\nx\t1\tysinglecopy\t1\t100\t1\t0\t100.0\n").unwrap();
        std::fs::write(&m, "#y_dose_scale=0.5\nclass\tK\ttech\tlength\na\t100\t\ta\n").unwrap();
        let mut o = vec![];
        super::estimate(s.to_str().unwrap(), m.to_str().unwrap(), "illumina", &mut o);
        std::fs::remove_dir_all(&d).unwrap();
        assert_eq!(String::from_utf8(o).unwrap(), "sample\tclass\testimate_bp\trelative\tmarker_sum_gc\tysinglecopy_sum\ty_dose\ty_status\nx\ta\t10\t0.1\t10.0\t100\t2.0000\ty_found\n");
    }
    #[test]
    fn gc_curve() {
        // pooled bin means recover a G+C-dependent loss at 0.5 counts per marker
        let mut x = 11u64;
        let mut by: Vec<Vec<u64>> = vec![vec![]; 32];
        let mut y = vec![];
        for i in 0..3_000_000u64 {
            let b = 10 + (i % 4) as usize;
            let p = 0.25 * (1.0 + 0.1 * (b as f64 - 11.5));
            let mut k = 0u64;
            for _ in 0..2 {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if ((x >> 11) as f64 / (1u64 << 53) as f64) < p {
                    k += 1;
                }
            }
            by[b].push(k);
            y.push(k);
        }
        let d = super::gc_bins(&y, &by);
        assert_eq!(d.len(), 4);
        for (b, r) in d {
            let e = 1.0 + 0.1 * (b as f64 - 11.5);
            assert!((r - e).abs() < 0.01, "bin {b}: {r} vs {e}");
        }
        assert!(super::gc_bins(&y[..1000], &vec![vec![]; 32]).is_empty()); // too few ruler counts
    }
    #[test]
    fn python_repr() {
        // values from python3 repr()
        for (x, s) in [
            (4.517420114170894e-05, "4.517420114170894e-05"),
            (0.0001756600924138287, "0.0001756600924138287"),
            (1.0095530193019127, "1.0095530193019127"),
            (27.94310690421051, "27.94310690421051"),
            (14.0, "14.0"),
            (1e16, "1e+16"),
            (123456789012345678.0, "1.2345678901234568e+17"),
            (0.0001, "0.0001"),
            (2.5e-10, "2.5e-10"),
            (14.5, "14.5"),
        ] {
            assert_eq!(super::py_float(x), s);
        }
    }
}

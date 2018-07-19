// tslint:disable:no-bitwise
export const makeArray = (length: number, val?: number | string) => {
  if (!val) {
    val = 0 | 0;
  }
  if (val < 0) {
    val = 0;
  }
  length |= 0;
  let max = 0;
  for (let i = length; i !== 0; i >>>= 1) {
    max++;
  }
  const n = Array(max);
  n[0] = [val];
  for (let i = 1; i < max; i++) {
    n[i] = n[i - 1].concat(n[i - 1]);
  }
  let a = new Array();
  for (let i = 0, l = length; l !== 0; l >>>= 1, i++) {
    if (l & 1) {
      a = a.concat(n[i]);
    }
  }
  return a;
};

export const bitCount = (() => {
  const a = new Uint8Array(256);
  let bin;
  for (let i = 0; i < 256; i++) {
    bin = i;
    bin = bin - ((bin >> 1) & 0x55555555);
    bin = (bin & 0x33333333) + ((bin >> 2) & 0x33333333);
    a[i] = (((bin + (bin >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
  }
  return a;
})();

export const nucleotideTo4Bit: { [key: string]: number } = {
  A: 8, // 0b1000
  C: 1, // 0b0001
  G: 2, // 0b0010
  T: 4, // 0b0100
};

export const setNucleotide = (...nucleotides: string[]) => {
  const n = nucleotides[0];
  nucleotideTo4Bit[n] = 0;
  for (let i = 1; i < nucleotides.length; i++) {
    nucleotideTo4Bit[n] |= nucleotideTo4Bit[nucleotides[i]];
  }
};

setNucleotide('-');
setNucleotide('W', 'A', 'T');
setNucleotide('S', 'G', 'C');
setNucleotide('M', 'A', 'C');
setNucleotide('K', 'G', 'T');
setNucleotide('R', 'A', 'G');
setNucleotide('Y', 'C', 'T');
setNucleotide('B', 'C', 'G', 'T');
setNucleotide('D', 'A', 'G', 'T');
setNucleotide('H', 'A', 'C', 'T');
setNucleotide('V', 'A', 'C', 'G');
setNucleotide('N', 'A', 'T', 'G', 'C');

// tslint:disable-next-line:variable-name
export const __4BitToNucleotide = (() => {
  const a = makeArray(16);
  const keys = Object.keys(nucleotideTo4Bit);
  for (let i = 0, len = keys.length; i < len; i++) {
    a[nucleotideTo4Bit[keys[i]]] = keys[i];
  }
  return a;
})();

export const nucleotideList = Object.keys(nucleotideTo4Bit);

// tslint:disable-next-line:variable-name
export const complementNucleotideHelper = (() => {
  const a = Object.create(null);
  a.A = 'T';
  a.G = 'C';
  a.B = 'V';
  a.H = 'D';
  a.M = 'K';
  a.R = 'Y';
  // S, W, N, - not included
  const keys = Object.keys(a);
  for (let i = 0, len = keys.length; i < len; i++) {
    a[a[keys[i]]] = keys[i];
  }
  a.S = 'S';
  a.W = 'W';
  a.N = 'N';
  a['-'] = '-';
  return a;
})();

export const nucleotideToBin = (s: string) => nucleotideTo4Bit[s] | 0;

export const binToNucleotide = (b: number) => __4BitToNucleotide[b] || '-';

export const complementNucleotide = (s: string) => complementNucleotideHelper[s];

export const NT = (() => {
  const a = Object.create(null);
  const b = new Uint8Array(256);
  const c = makeArray(256);
  const d = makeArray(256);

  const keys = Object.keys(nucleotideTo4Bit);
  const len = keys.length;
  let ki;
  let kj;
  let byte;

  for (let i = 0; i < len; i++) {
    ki = keys[i];
    for (let j = 0; j < len; j++) {
      kj = keys[j];
      byte = nucleotideTo4Bit[ki] | (nucleotideTo4Bit[kj] << 4);
      a[ki + kj] = byte;
      b[byte] = nucleotideTo4Bit[complementNucleotide(kj)] | (nucleotideTo4Bit[complementNucleotide(ki)] << 4);
      c[byte] = ki + kj;
      d[byte] = [ki, kj];
    }
  }

  return {
    byteComplement: b,
    byteNucleotideContent: d,
    byteToNucleotides: c,
    nucleotidesToByteHelper: a,
  };
})();

export const nucleotidesToByte = (ss: string) => NT.nucleotidesToByteHelper[ss] | 0;

/* amino acids */

export const codonToAminoAcid: { [key: string]: string } = {
  AAA: 'K',
  AAC: 'N',
  AAG: 'K',
  AAT: 'N',
  ACA: 'T',
  ACC: 'T',
  ACG: 'T',
  ACT: 'T',
  AGA: 'R',
  AGC: 'S',
  AGG: 'R',
  AGT: 'S',
  ATA: 'I',
  ATC: 'I',
  ATG: 'M',
  ATT: 'I',
  CAA: 'Q',
  CAC: 'H',
  CAG: 'Q',
  CAT: 'H',
  CCA: 'P',
  CCC: 'P',
  CCG: 'P',
  CCT: 'P',
  CGA: 'R',
  CGC: 'R',
  CGG: 'R',
  CGT: 'R',
  CTA: 'L',
  CTC: 'L',
  CTG: 'L',
  CTT: 'L',
  GAA: 'E',
  GAC: 'D',
  GAG: 'E',
  GAT: 'D',
  GCA: 'A',
  GCC: 'A',
  GCG: 'A',
  GCT: 'A',
  GGA: 'G',
  GGC: 'G',
  GGG: 'G',
  GGT: 'G',
  GTA: 'V',
  GTC: 'V',
  GTG: 'V',
  GTT: 'V',
  TAA: '*',
  TAC: 'Y',
  TAG: '&',
  TAT: 'Y',
  TCA: 'S',
  TCC: 'S',
  TCG: 'S',
  TCT: 'S',
  TGA: '$',
  TGC: 'C',
  TGG: 'W',
  TGT: 'C',
  TTA: 'L',
  TTC: 'F',
  TTG: 'L',
  TTT: 'F',
};

// tslint:disable-next-line:variable-name
export const __12BitToAminoAcid = (() => {
  const a = makeArray(4096, '?');
  const codons = Object.keys(codonToAminoAcid);
  let codon;
  for (let i = 0, len = codons.length; i < len; i++) {
    codon = codons[i];
    a[(nucleotideToBin(codon[2]) << 8) | (nucleotideToBin(codon[1]) << 4) | nucleotideToBin(codon[0])] =
      codonToAminoAcid[codon];
  }

  return a;
})();

export default NT;

"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
// tslint:disable:no-bitwise
exports.makeArray = function (length, val) {
    if (!val) {
        val = 0 | 0;
    }
    if (val < 0) {
        val = 0;
    }
    length |= 0;
    var max = 0;
    for (var i = length; i !== 0; i >>>= 1) {
        max++;
    }
    var n = Array(max);
    n[0] = [val];
    for (var i = 1; i < max; i++) {
        n[i] = n[i - 1].concat(n[i - 1]);
    }
    var a = new Array();
    for (var i = 0, l = length; l !== 0; l >>>= 1, i++) {
        if (l & 1) {
            a = a.concat(n[i]);
        }
    }
    return a;
};
exports.makeBitCount = function () {
    var a = new Uint8Array(256);
    var bin;
    for (var i = 0; i < 256; i++) {
        bin = i;
        bin = bin - ((bin >> 1) & 0x55555555);
        bin = (bin & 0x33333333) + ((bin >> 2) & 0x33333333);
        a[i] = (((bin + (bin >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
    }
    return a;
};
exports.nucleotideTo4Bit = {
    A: 8,
    C: 1,
    G: 2,
    T: 4,
};
exports.setNucleotide = function () {
    var nucleotides = [];
    for (var _i = 0; _i < arguments.length; _i++) {
        nucleotides[_i] = arguments[_i];
    }
    var n = nucleotides[0];
    exports.nucleotideTo4Bit[n] = 0;
    for (var i = 1; i < nucleotides.length; i++) {
        exports.nucleotideTo4Bit[n] |= exports.nucleotideTo4Bit[nucleotides[i]];
    }
};
exports.setNucleotide('-');
exports.setNucleotide('W', 'A', 'T');
exports.setNucleotide('S', 'G', 'C');
exports.setNucleotide('M', 'A', 'C');
exports.setNucleotide('K', 'G', 'T');
exports.setNucleotide('R', 'A', 'G');
exports.setNucleotide('Y', 'C', 'T');
exports.setNucleotide('B', 'C', 'G', 'T');
exports.setNucleotide('D', 'A', 'G', 'T');
exports.setNucleotide('H', 'A', 'C', 'T');
exports.setNucleotide('V', 'A', 'C', 'G');
exports.setNucleotide('N', 'A', 'T', 'G', 'C');
// tslint:disable-next-line:variable-name
exports.__4BitToNucleotide = (function () {
    var a = exports.makeArray(16);
    var keys = Object.keys(exports.nucleotideTo4Bit);
    for (var i = 0, len = keys.length; i < len; i++) {
        a[exports.nucleotideTo4Bit[keys[i]]] = keys[i];
    }
    return a;
})();
exports.nucleotideList = Object.keys(exports.nucleotideTo4Bit);
// tslint:disable-next-line:variable-name
exports.complementNucleotideHelper = (function () {
    var a = Object.create(null);
    a.A = 'T';
    a.G = 'C';
    a.B = 'V';
    a.H = 'D';
    a.M = 'K';
    a.R = 'Y';
    // S, W, N, - not included
    var keys = Object.keys(a);
    for (var i = 0, len = keys.length; i < len; i++) {
        a[a[keys[i]]] = keys[i];
    }
    a.S = 'S';
    a.W = 'W';
    a.N = 'N';
    a['-'] = '-';
    return a;
})();
exports.nucleotideToBin = function (s) { return exports.nucleotideTo4Bit[s] | 0; };
exports.binToNucleotide = function (b) { return exports.__4BitToNucleotide[b] || '-'; };
exports.complementNucleotide = function (s) { return exports.complementNucleotideHelper[s]; };
exports.NT = (function () {
    var a = Object.create(null);
    var b = new Uint8Array(256);
    var c = exports.makeArray(256);
    var d = exports.makeArray(256);
    var keys = Object.keys(exports.nucleotideTo4Bit);
    var len = keys.length;
    var ki;
    var kj;
    var byte;
    for (var i = 0; i < len; i++) {
        ki = keys[i];
        for (var j = 0; j < len; j++) {
            kj = keys[j];
            byte = exports.nucleotideTo4Bit[ki] | (exports.nucleotideTo4Bit[kj] << 4);
            a[ki + kj] = byte;
            b[byte] = exports.nucleotideTo4Bit[exports.complementNucleotide(kj)] | (exports.nucleotideTo4Bit[exports.complementNucleotide(ki)] << 4);
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
exports.nucleotidesToByte = function (ss) { return exports.NT.nucleotidesToByteHelper[ss] | 0; };
/* amino acids */
exports.codonToAminoAcid = {
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
exports.__12BitToAminoAcid = (function () {
    var a = exports.makeArray(4096, '?');
    var codons = Object.keys(exports.codonToAminoAcid);
    var codon;
    for (var i = 0, len = codons.length; i < len; i++) {
        codon = codons[i];
        a[(exports.nucleotideToBin(codon[2]) << 8) | (exports.nucleotideToBin(codon[1]) << 4) | exports.nucleotideToBin(codon[0])] =
            exports.codonToAminoAcid[codon];
    }
    return a;
})();
exports.default = exports.NT;
//# sourceMappingURL=nt.js.map
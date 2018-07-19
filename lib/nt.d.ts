export declare const makeArray: (length: number, val?: string | number | undefined) => any[];
export declare const makeBitCount: () => Uint8Array;
export declare const nucleotideTo4Bit: {
    [key: string]: number;
};
export declare const setNucleotide: (...nucleotides: string[]) => void;
export declare const __4BitToNucleotide: any[];
export declare const nucleotideList: string[];
export declare const complementNucleotideHelper: any;
export declare const nucleotideToBin: (s: string) => number;
export declare const binToNucleotide: (b: number) => any;
export declare const complementNucleotide: (s: string) => any;
export declare const NT: {
    byteComplement: Uint8Array;
    byteNucleotideContent: any[];
    byteToNucleotides: any[];
    nucleotidesToByteHelper: any;
};
export declare const nucleotidesToByte: (ss: string) => number;
export declare const codonToAminoAcid: {
    [key: string]: string;
};
export declare const __12BitToAminoAcid: any[];
export default NT;

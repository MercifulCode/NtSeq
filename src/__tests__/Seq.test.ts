import { Seq } from '../Seq';
import * as Path from 'path';

describe('Seq', () => {
  it('Should report the correct sequence size.', () => {
    const seq = new Seq().read('ATGCATG');
    expect(seq.size).toEqual(7);
  });

  it('Should report the correct sequence.', () => {
    const expected = 'ATGCATG';
    const seq = new Seq().read(expected);
    expect(seq.sequence()).toEqual(expected);
  });

  it('Should report the correct sequence compliment.', () => {
    const expected = 'CATGCAT';
    const seq = new Seq().read('ATGCATG');
    expect(seq.complement.sequence()).toEqual(expected);
  });

  it('Should report sequence equivalence.', () => {
    const seqA = new Seq().read('ATGCATGC');
    const seqB = new Seq().read('ATGCATGC');
    expect(seqA.equivalent(seqB)).toBeTruthy();
  });

  it('Should report sequence equivalence.', () => {
    const seqA = new Seq().read('ATGCATGC');
    const seqB = new Seq().read('TAGCATGC');
    expect(seqA.equivalent(seqB)).not.toBeTruthy();
  });

  it('Should replicate the sequence correctly.', () => {
    const expected = 'ATGCATGCA';
    const seqA = new Seq().read(expected);
    expect(seqA.replicate().sequence()).toEqual(expected);
  });

  it('Should replicate with one parameter (offset 1) for a sequence.', () => {
    const seq = new Seq().read('ATGCATGCA');
    const expected = 'ATGCATGC';
    expect(seq.replicate(0, 8).sequence()).toEqual(expected);
  });

  it('Should replicate with subset of sequence.', () => {
    const seq = new Seq().read('TGCGCTACGAAGACGT');
    const expected = 'GCGCTACG';
    expect(seq.replicate(1, 8).sequence()).toEqual(expected);
  });

  it('Should perform polymerization correctly.', () => {
    const seqA = new Seq().read('ATGC');
    const seqB = new Seq().read('TCAG');
    const expected = 'ATGCTCAG';
    expect(seqA.polymerize(seqB).sequence()).toEqual(expected);
  });

  it('Should perform sequence insertion correctly.', () => {
    const seqA = new Seq().read('ATGC');
    const seqB = new Seq().read('TCAG');
    const expected = 'TCAGATGC';
    expect(seqA.insertion(seqB, 0).sequence()).toEqual(expected);
  });

  it('Should perform sequence insertion correctly when given an offset.', () => {
    const seqA = new Seq().read('ATGC');
    const seqB = new Seq().read('TCAG');
    const expected = 'ATCAGTGC';
    expect(seqA.insertion(seqB, 1).sequence()).toEqual(expected);
  });

  it('Should delete when given a offset.', () => {
    const seq = new Seq().read('ATGCATGCA');
    const expected = 'ACATGCA';
    expect(seq.deletion(1, 2).sequence()).toEqual(expected);
  });

  it('Should delete when not given a offset.', () => {
    const seq = new Seq().read('ATGCATGCA');
    const expected = 'TGCATGCA';
    expect(seq.deletion(0, 1).sequence()).toEqual(expected);
  });

  it('Should allow a sequence to be repeated.', () => {
    const seq = new Seq().read('ATG');
    const expected = 'ATGATG';
    expect(seq.repeat(2).sequence()).toEqual(expected);
  });

  it('Should allow a sequence to be repeated a lot.', () => {
    const seq = new Seq().read('ATG');
    const expected = 'ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG';
    expect(seq.repeat(17).sequence()).toEqual(expected);
  });

  it('Should allow a sequence to be masked.', () => {
    const seqA = new Seq().read('ATGCATGC');
    const seqB = new Seq().read('AACCATNC');
    const expected = 'A--CATGC';
    expect(seqA.mask(seqB).sequence()).toEqual(expected);
  });

  it('Should allow a sequence to be covered.', () => {
    const seqA = new Seq().read('ATGCATGC');
    const seqB = new Seq().read('AACCATNC');
    const expected = 'AWSCATNC';
    expect(seqA.cover(seqB).sequence()).toEqual(expected);
  });

  it('Should have the correct content for a sequence.', () => {
    const seq = new Seq().read('ATGCATGCATGC');
    const { content } = seq;
    expect(content['A']).toEqual(3);
    expect(content['C']).toEqual(3);
    expect(content['G']).toEqual(3);
    expect(content['T']).toEqual(3);
  });

  it('Should have the correct content for a sequence with degenerate nucleotides.', () => {
    const seq = new Seq().read('AT-GCATGCATGCATGC--NW');
    const { content } = seq;
    expect(content['A']).toEqual(4);
    expect(content['C']).toEqual(4);
    expect(content['G']).toEqual(4);
    expect(content['T']).toEqual(4);
    expect(content['-']).toEqual(3);
    expect(content['N']).toEqual(1);
    expect(content['W']).toEqual(1);
  });

  it('Should have the correct fractional content for a sequence.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    
    const { fractionalContent } = seq;
    expect(fractionalContent['A']).toEqual(0.25);
    expect(fractionalContent['C']).toEqual(0.25);
    expect(fractionalContent['G']).toEqual(0.25);
    expect(fractionalContent['T']).toEqual(0.25);
  });

  it('Should have the correct ATGC (average) content for a sequence.', () => {
    const seq = new Seq().read('WWWWNNNN');
    
    const { contentATGC } = seq;
    expect(contentATGC['A']).toEqual(3);
    expect(contentATGC['C']).toEqual(1);
    expect(contentATGC['G']).toEqual(1);
    expect(contentATGC['T']).toEqual(3);
  });

  it('Should have the correct ATGC (average) fractional content for a sequence.', () => {
    const seq = new Seq().read('WWWWNNNN');
    
    const { fractionalContentATGC } = seq;
    expect(fractionalContentATGC['A']).toEqual(0.375);
    expect(fractionalContentATGC['C']).toEqual(0.125);
    expect(fractionalContentATGC['G']).toEqual(0.125);
    expect(fractionalContentATGC['T']).toEqual(0.375);
  });

  it('Should translate a sequence consisting of a single amino acid.', () => {
    const seq = new Seq().read('ATG');
    const expected = 'M';
    expect(seq.translate()).toEqual(expected);
  });

  it('Should translate a sequence with several amino acids.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'MHACM';
    expect(seq.translate()).toEqual(expected);
  });

  it('Should translate a sequence with several amino acids and a given offset.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'CMHAC';
    expect(seq.translate(1)).toEqual(expected);
  });

  it('Should translate a sequence with several amino acids and both a given offset and length.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'CM';
    expect(seq.translate(1, 7)).toEqual(expected);
  });

  it('Should translate a frame correctly.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'MHACM';
    expect(seq.translateFrame()).toEqual(expected);
  });

  it('Should translate a frame correctly when given an offset.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'HACM';
    expect(seq.translateFrame(0, 1)).toEqual(expected);
  });

  it('Should translate a frame correctly when given an offset and length.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'HA';
    expect(seq.translateFrame(0, 1, 2)).toEqual(expected);
  });

  it('Should translate a frame correctly when given a specific frame, offset and length.', () => {
    const seq = new Seq().read('ATGCATGCATGCATGC');
    const expected = 'MH';
    expect(seq.translateFrame(1, 1, 2)).toEqual(expected);
  });

  it.each(['fasta', '4bnt'])(
    'Should have mapSequence successfully find match scores loaded from .%s file.', (ext) => {
    const seqA = new Seq().loadFile(Path.resolve(__dirname, `data/sequence.${ext}`), ext)
    const seqB = new Seq().read('TCTTATTTGTGCTGTTTATT');
    const matchMap = seqA.mapSequence(seqB).initialize().sort();
    const testData = [
      233,
      1571,
      5942,
      14035,
      24226,
      31360,
      32021,
      26118,
      17738,
      9439,
      4181,
      1509,
      418,
      112,
      15,
      3,
      0,
      0,
      0,
      0,
      1
    ];

    const matchFrequencyData = matchMap.matchFrequencyData();
    for (let i = 0; i < matchFrequencyData.length; i++) {
      expect(matchFrequencyData[i]).toEqual(testData[i]);
    }
    expect(matchMap.best().position).toEqual(34786);
    expect(matchMap.best().score).toEqual(20);
  });
})
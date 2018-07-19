// tslint:disable:no-bitwise
import { MatchMap } from './MatchMap';
import * as nt from './nt';

export class Seq {
  public buffer = new ArrayBuffer(4);
  public endPadding = 0;
  public length = 0;

  protected myComplement = null as null | ArrayBuffer;
  protected myContent = null as any;
  protected myContentATGC = null as any;
  protected myFractionalContent = null as any;
  protected myFractionalContentATGC?: any;

  public get isRNA() {
    return this.type === 'RNA';
  }

  public get size() {
    return this.length;
  }

  constructor(readonly type: string = 'DNA') {
    // TODO Wat
    /*
    if (!{ RNA: true, DNA: true }[type]) {
      throw new Error('Sequence type ' + type + ' not supported');
    }
    */
    if (type !== 'DNA' && type !== 'RNA') {
      throw new Error('Sequence type ' + type + ' not supported');
    }
  }

  public read(strData: string) {
    const ntToByte = nt.nucleotidesToByte;

    const nucleotideString = strData
      .toUpperCase()
      .replace(/\s/g, '')
      .replace('U', 'T')
      .replace(/[^ATGCBVHDMKRYSWN\-]/g, '-');

    const length = nucleotideString.length | 0;
    const max = length >>> 1;
    const odd = length & 1;

    const endPadding = (4 - ((max + odd) % 4)) % 4;
    this.endPadding = endPadding;

    const buffer = new ArrayBuffer(max + odd + endPadding + 4);
    const dataArray = new Int8Array(buffer, 4);

    let n;
    let i = 0;

    for (; i < max; i++) {
      n = i << 1;
      dataArray[i] = ntToByte(nucleotideString[n] + nucleotideString[++n]);
    }

    if (odd) {
      dataArray[i] = nt.nucleotideTo4Bit[nucleotideString[i << 1]];
    }

    this.buffer = buffer;
    this.length = length;
    new Uint32Array(buffer, 0, 1)[0] = length;

    this.myComplement = null;

    this.myContent = null;
    this.myFractionalContent = null;
    this.myContentATGC = null;
    this.myFractionalContent = null;

    return this;
  }

  public readBuffer(buffer: ArrayBuffer) {
    this.buffer = buffer;

    const length = new Uint32Array(buffer, 0, 1)[0];

    const max = length >>> 1;
    const odd = length & 1;

    const endPadding = (4 - ((max + odd) % 4)) % 4;

    this.endPadding = endPadding;
    this.length = length;
    this.myComplement = null;
    this.myContent = null;
    this.myFractionalContent = null;
    this.myContentATGC = null;
    this.myFractionalContent = null;

    return this;
  }

  public readFASTA(strFASTA: string) {
    const data = strFASTA.split(/\n\r?/gi);

    while (data.length && data[0][0] === '>') {
      data.shift();
    }

    return this.read(data.join(''));
  }

  public byteComplement() {
    const bComp = nt.NT.byteComplement;

    const fwdBuffer = this.buffer;
    const len = fwdBuffer.byteLength;

    let n;
    let i;

    let copyBuffer;
    let fromArray;
    let copyArray;

    const isOdd = this.length & 1;

    if (isOdd) {
      copyBuffer = new ArrayBuffer(len);
      copyArray = new Uint32Array(copyBuffer, 4);
      fromArray = new Uint32Array(fwdBuffer, 4);

      n = (len - 4) >>> 2;
      while (n--) {
        copyArray[n] = (fromArray[n] << 4) | (fromArray[n - 1] >>> 28);
      }
    } else {
      copyBuffer = fwdBuffer;
    }

    const fwdArray = new Uint8Array(copyBuffer, 4);

    const buffer = new ArrayBuffer(len);
    const dataArray = new Uint8Array(buffer, 4);

    n = len - 4 - this.endPadding;
    i = 0;
    while (n--) {
      dataArray[i++] = bComp[fwdArray[n]];
    }

    new Uint32Array(buffer, 0, 1)[0] = this.length;

    return buffer;
  }

  public sequence() {
    const byteToNt = nt.NT.byteToNucleotides;
    const buffer = this.buffer;

    if (buffer.byteLength < 4) {
      return '';
    }

    const dataArray = new Uint8Array(buffer, 4);
    const len = buffer.byteLength - 4 - this.endPadding;

    const nts = nt.makeArray(len);

    let i = 0;
    do {
      nts[i] = byteToNt[dataArray[i]];
      ++i;
    } while (i < len);

    let returnString;

    i = nts.length - 1;
    if (this.length & 1) {
      nts[i] = nts[i][0];
    }

    returnString = nts.join('');

    if (this.isRNA) {
      returnString = returnString.replace(/T/gi, 'U');
    }

    return returnString;
  }

  public get complement() {
    if (!this.myComplement) {
      this.myComplement = this.byteComplement();
    }

    const complement = new Seq(this.type).readBuffer(this.myComplement.slice(0));

    complement.myComplement = this.buffer.slice(0);

    return complement;
  }

  public equivalent(seq: Seq) {
    if (!(seq instanceof Seq)) {
      throw new Error('Can only check for equivalence between sequences');
    }

    if (this.type !== seq.type) {
      return false;
    }

    const checkInts = new Uint32Array(this.buffer);
    const compareInts = new Uint32Array(seq.buffer);

    for (let i = 0, len = checkInts.length; i < len; i++) {
      if (checkInts[i] !== compareInts[i]) {
        return false;
      }
    }

    return true;
  }

  public replicate(start: number = 0, length?: number) {
    start |= 0;

    if (start < 0) {
      start = this.length + start;
    }

    if (length === undefined) {
      if (start === 0) {
        return this.clone();
      }

      length = this.length - start;
    } else {
      length |= 0;
      length = Math.min(length, this.length - start);
    }

    length = Math.min(length, this.length - start);

    if (length <= 0) {
      return this.nullSeq();
    }

    return this.slice(start, length);
  }

  public polymerize(seq: Seq) {
    const seqLen = seq.length;

    if (!(seq instanceof Seq)) {
      throw new Error('.polymerize requires valid sequence');
    }

    if (!this.length) {
      return seq.clone();
    }

    const length = this.length + seqLen;

    const max = length >>> 1;
    const odd = length & 1;

    const endPadding = (4 - ((max + odd) % 4)) % 4;
    const newBuffer = new ArrayBuffer(max + odd + endPadding + 4);
    const newArray = new Uint32Array(newBuffer, 4);

    const copyBuffer = this.buffer;
    const copyArray = new Uint32Array(copyBuffer, 4);

    const seqBuffer = seq.buffer;
    const seqArray = new Uint32Array(seqBuffer, 4);

    let copyPos = 0;
    const shift = (this.length % 8) * 4;
    const shiftSeq = 32 - shift;

    for (const len = copyArray.length; copyPos < len; copyPos++) {
      newArray[copyPos] = copyArray[copyPos];
    }

    if (shift) {
      newArray[--copyPos] |= seqArray[0] << shift;

      for (let i = 0, len = seqArray.length; i < len; i++) {
        newArray[++copyPos] = (seqArray[i] >>> shiftSeq) | (seqArray[i + 1] << shift);
      }
    } else {
      for (let i = 0, len = seqArray.length; i < len; i++) {
        newArray[copyPos++] = seqArray[i];
      }
    }

    new Uint32Array(newBuffer, 0, 1)[0] = length;

    return new Seq(this.type).readBuffer(newBuffer);
  }

  public insertion(seq: Seq, offset: number) {
    if (!(seq instanceof Seq)) {
      throw new Error('Insertion requires valid sequence');
    }

    offset |= 0;

    if (offset < 0) {
      offset = this.length + offset;
    }

    offset = Math.min(offset, this.length);

    return this.replicate(0, offset)
      .polymerize(seq)
      .polymerize(this.replicate(offset));
  }

  public deletion(offset: number, count: number) {
    if (offset === undefined || count === undefined) {
      throw new Error('Must give valid offset and count for deletion');
    }

    offset |= 0;
    count |= 0;

    if (count === 0) {
      return this.clone();
    }

    if (count < 0) {
      throw new Error('Invalid count for deletion');
    }

    if (offset < 0) {
      offset = this.length + offset;
    }

    offset = Math.min(offset, this.length);

    return this.replicate(0, offset).polymerize(this.replicate(offset + count));
  }

  public repeat(count: number) {
    count |= 0;

    let copy = this.replicate();
    let base = new Seq(this.type);

    if (count <= 0) {
      return base;
    }

    while (true) {
      if (count & 1) {
        base = base.polymerize(copy);
      }
      count >>>= 1;
      if (!count) {
        break;
      }
      copy = copy.polymerize(copy);
    }

    return base;
  }

  public mask(seq: Seq) {
    if (!(seq instanceof Seq)) {
      throw new Error('Can only mask with valid sequence');
    }

    const newBuffer = this.buffer.slice(0);
    const newArray = new Uint32Array(newBuffer, 4);
    const compareArray = new Uint32Array(seq.buffer, 4);

    for (let i = 0, len = newArray.length; i < len; i++) {
      newArray[i] &= compareArray[i];
    }

    return new Seq(this.type).readBuffer(newBuffer);
  }

  public cover(seq: Seq) {
    if (!(seq instanceof Seq)) {
      throw new Error('Can only cover with valid sequence');
    }

    const newBuffer = this.buffer.slice(0);
    const newArray = new Uint32Array(newBuffer, 4);
    const compareArray = new Uint32Array(seq.buffer, 4);

    for (let i = 0, len = newArray.length; i < len; i++) {
      newArray[i] |= compareArray[i];
    }

    return new Seq(this.type).readBuffer(newBuffer);
  }

  public nullSeq() {
    return new Seq(this.type).readBuffer(new ArrayBuffer(4));
  }

  public clone() {
    return new Seq(this.type).readBuffer(this.buffer.slice(0));
  }

  public slice(start: number, length: number) {
    const max = length >>> 1;
    const odd = length & 1;

    const endPadding = (4 - ((max + odd) % 4)) % 4;
    const newBuffer = new ArrayBuffer(max + odd + endPadding + 4);
    const newArray = new Uint32Array(newBuffer, 4);

    const subBuffer = this.buffer.slice(4 + (start >>> 1), 4 + (start >>> 1) + newBuffer.byteLength);
    const subInt32Length = subBuffer.byteLength >>> 2;
    const subArray = new Uint32Array(subBuffer, 0, subInt32Length);

    if (start & 1) {
      let subArrayIndex = 0;
      do {
        newArray[subArrayIndex] = (subArray[subArrayIndex] >>> 4) | (subArray[subArrayIndex + 1] << 28);
        ++subArrayIndex;
      } while (subArrayIndex < subArray.length);

      const remainder = subBuffer.byteLength - subArray.byteLength;
      if (remainder) {
        const remainderArray = new Uint8Array(newBuffer, 4 + (subArrayIndex << 2));
        const subRemainderArray = new Uint8Array(subBuffer, subArrayIndex << 2);
        if (newArray.length > 0) {
          newArray[subArrayIndex - 1] |= subRemainderArray[0] << 28;
        }
        for (let i = 0, len = subRemainderArray.length; i < len; i++) {
          remainderArray[i] = (subRemainderArray[i] >>> 4) | (subRemainderArray[i + 1] << 4);
        }
      }
    } else {
      let subArrayIndex = 0;
      do {
        newArray[subArrayIndex] = subArray[subArrayIndex];
        ++subArrayIndex;
      } while (subArrayIndex < subArray.length);

      const remainder = subArray.byteLength - subBuffer.byteLength;
      if (remainder) {
        const remainderArray = new Uint8Array(newBuffer, 4 + (subArrayIndex << 2));
        const subRemainderArray = new Uint8Array(subBuffer, subArrayIndex << 2);
        for (let i = 0, len = subRemainderArray.length; i < len; i++) {
          remainderArray[i] = subRemainderArray[i];
        }
      }
    }

    const clearShift = (endPadding * 2 + odd) * 4;

    const clearOut = new Uint32Array(newBuffer, newBuffer.byteLength - 4);
    clearOut[0] = (clearOut[0] << clearShift) >>> clearShift;

    new Uint32Array(newBuffer, 0, 1)[0] = length;

    return new Seq(this.type).readBuffer(newBuffer);
  }

  public get content() {
    if (!this.content) {
      const ntContentByte = nt.makeArray(256);

      const buffer = this.buffer;
      const dataArray = new Uint8Array(buffer);

      for (let i = 4; i < buffer.byteLength - this.endPadding; i++) {
        ntContentByte[dataArray[i]]++;
      }

      const binToNt = nt.binToNucleotide;
      const ntList = nt.nucleotideList;

      const ntContent = Object.create(null);
      for (let i = 0, len = ntList.length; i < len; i++) {
        ntContent[ntList[i]] = 0;
      }

      for (let i = 0, len = ntContentByte.length; i < len; i++) {
        if (ntContentByte[i]) {
          ntContent[binToNt(i & 0xf)] += ntContentByte[i];
          ntContent[binToNt(i >>> 4)] += ntContentByte[i];
        }
      }

      if (this.length & 1) {
        ntContent['-']--;
      }

      if (this.isRNA) {
        ntContent.U = ntContent.T;
        delete ntContent.T;
      }

      this.myContent = ntContent;
    }

    const returnContent = Object.create(null);
    const keys = Object.keys(this.content);
    for (let i = 0, len = keys.length; i < len; i++) {
      returnContent[keys[i]] = this.content[keys[i]];
    }

    return returnContent;
  }

  public get fractionalContent() {
    if (!this.myFractionalContent) {
      const content = this.content();
      const nts = Object.keys(content);
      for (let i = 0, len = nts.length; i < len; i++) {
        content[nts[i]] = content[nts[i]] / this.length;
      }

      this.myFractionalContent = content;
    }

    const returnContent = Object.create(null);
    const keys = Object.keys(this.myFractionalContent);
    for (let i = 0, len = keys.length; i < len; i++) {
      returnContent[keys[i]] = this.myFractionalContent[keys[i]];
    }

    return returnContent;
  }

  public get contentATGC() {
    if (!this.myContentATGC) {
      const ntToBin = nt.nucleotideToBin;

      const content = this.content();
      const nts = Object.keys(content);
      const contentATGC = Object.create(null);
      contentATGC.A = 0;
      contentATGC.T = 0;
      contentATGC.G = 0;
      contentATGC.C = 0;

      let bits = 0;
      let nucleotides;
      let ntBin;
      let n;
      let curContent;

      for (let i = 0, len = nts.length; i < len; i++) {
        nucleotides = nts[i];
        n = ntToBin(nucleotides);
        for (bits = 0; n; bits++) {
          n &= n - 1;
        }

        ntBin = ntToBin(nucleotides);
        curContent = content[nts[i]] * (1 / bits);

        contentATGC.A += (ntToBin('A') & ntBin) | 0 && curContent;
        contentATGC.T += (ntToBin('T') & ntBin) | 0 && curContent;
        contentATGC.G += (ntToBin('G') & ntBin) | 0 && curContent;
        contentATGC.C += (ntToBin('C') & ntBin) | 0 && curContent;
      }

      if (this.isRNA) {
        contentATGC.U = contentATGC.T;
        delete contentATGC.T;
      }

      this.myContentATGC = contentATGC;
    }

    const returnContent = Object.create(null);
    const keys = Object.keys(this.myContentATGC);
    for (let i = 0, len = keys.length; i < len; i++) {
      returnContent[keys[i]] = this.myContentATGC[keys[i]];
    }

    return returnContent;
  }

  public get fractionalContentATGC() {
    if (!this.fractionalContentATGC) {
      const content = this.myContentATGC();
      const nts = Object.keys(content);
      for (let i = 0, len = nts.length; i < len; i++) {
        content[nts[i]] = content[nts[i]] / this.length;
      }

      this.myFractionalContentATGC = content;
    }

    const returnContent = Object.create(null);
    const keys = Object.keys(this.fractionalContentATGC);
    for (let i = 0, len = keys.length; i < len; i++) {
      returnContent[keys[i]] = this.myFractionalContentATGC[keys[i]];
    }

    return returnContent;
  }

  public translate(ntOffset: number = 0, ntCount?: number) {
    const binToAA = nt.__12BitToAminoAcid;

    ntOffset |= 0;

    if (ntCount === undefined) {
      ntCount = this.length - ntOffset;
    }

    ntCount |= 0;
    ntCount -= ntCount % 3;

    const offset = (ntOffset >>> 1) + 4;
    const max = offset + (ntCount >>> 1) + (ntCount & 1);

    const dataArray = new Uint8Array(this.buffer);

    const aminoAcids = nt.makeArray(ntCount / 3);
    /**/
    let aa = 0;
    let lastByte;
    let byte1;
    let byte2;
    let byte3;

    if ((ntOffset & 1) === 0) {
      for (let i = offset; i < max; i += 3) {
        byte1 = dataArray[i];
        byte2 = dataArray[i + 1];
        byte3 = dataArray[i + 2];

        aminoAcids[aa++] = binToAA[byte1 | ((byte2 & 0xf) << 8)];
        aminoAcids[aa++] = binToAA[(byte3 << 4) | (byte2 >>> 4)];
      }
    } else {
      lastByte = dataArray[offset];

      for (let i = offset + 1; i < max; i += 3) {
        byte1 = dataArray[i];
        byte2 = dataArray[i + 1];
        byte3 = dataArray[i + 2];

        aminoAcids[aa++] = binToAA[(lastByte >> 4) | (byte1 << 4)];
        aminoAcids[aa++] = binToAA[byte2 | ((byte3 & 0xf) << 8)];

        lastByte = byte3;
      }
    }

    if (ntCount & 1) {
      aminoAcids.pop();
    }

    return aminoAcids.join('');
  }

  public translateFrame(frame: number = 0, AAoffset?: number, AAcount?: number) {
    if (frame !== 0 && frame !== 1 && frame !== 2) {
      throw new Error('Invalid translation frame, must be 0, 1 or 2.');
    }

    if (AAoffset === undefined) {
      return this.translate(frame);
    }

    if (AAcount === undefined) {
      return this.translate(frame + (AAoffset | 0) * 3);
    }

    return this.translate(frame + (AAoffset | 0) * 3, (AAcount * 3) | 0);
  }

  public mapSequence(seq: Seq, offset?: number) {
    if (!(seq instanceof Seq)) {
      throw new Error('.mapSequence requires valid Seq');
    }

    return new MatchMap(seq, this, offset);
  }
}

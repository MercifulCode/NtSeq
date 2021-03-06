import { MatchResult } from './MatchResult';
import { makeArray, makeBitCount } from './nt';
import { Seq } from './Seq';

// tslint:disable:no-bitwise
export class MatchMap {
  public myResults? = new Array<number>();
  public orderedResults = new Array<any>();
  public myMatchFrequencyData = null as null | any[];
  public initialized = false;
  public offset = 0;
  public positionAdjustment = 0;
  public debug = {
    prepareTime: null as null | number,
    searchTime: null as null | number,
    sortTime: null as null | number,
  };
  constructor(readonly query: Seq, readonly searchSpace: Seq, offset?: any) {}

  public initialize(results?: number[]) {
    this.orderedResults = [];
    this.myMatchFrequencyData = null;
    this.initialized = true;

    let t = new Date().valueOf();

    if (!results) {
      const dataArray = new Uint32Array(this.execute(this.query.buffer, this.searchSpace.buffer));

      this.debug.searchTime = -t + (t = new Date().valueOf());

      const queryLen = this.query.size;
      const searchLen = this.searchSpace.size;
      const resultsLen = ((8 - (queryLen % 8)) % 8) + 1;
      const totalLen = searchLen + queryLen - 1;
      results = [].slice.call(dataArray, resultsLen, resultsLen + totalLen);
    }

    this.myResults = results;
    this.debug.prepareTime = -t + (t = new Date().valueOf());

    return this;
  }

  public sort() {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    const t = new Date().valueOf();

    if (this.debug.sortTime !== null) {
      return this;
    }

    const adjust = this.positionAdjustment;
    this.orderedResults = this.myResults!.map((v, i) => ({ n: i + adjust, s: v })).sort((a, b) => b.s - a.s);
    this.debug.sortTime = new Date().valueOf() - t;
    return this;
  }

  public __calculate_p_match(query: Seq, searchSpace: Seq) {
    /*
    The approximate probability that two randomly chosen nucleotides
    from QUERY and SEARCH match each other
  */

    const queryContent = query.fractionalContentATGC();
    const searchSpaceContent = searchSpace.fractionalContentATGC();

    return (
      queryContent.A * searchSpaceContent.A +
      queryContent.T * searchSpaceContent.T +
      queryContent.G * searchSpaceContent.G +
      queryContent.C * searchSpaceContent.C
    );
  }

  public results(offset: number, count: number) {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    if (offset === undefined) {
      return this.myResults!.slice();
    }

    if (count === undefined) {
      return this.myResults!.slice(offset | 0);
    }

    return this.myResults!.slice(offset | 0, count | 0);
  }

  public best() {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    if (!this.orderedResults.length) {
      throw new Error('MatchMap must be sorted first.');
    }

    const result = this.orderedResults[0];
    return new MatchResult(this, result.n, result.s);
  }

  public top(n: number) {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    if (!this.orderedResults.length) {
      throw new Error('MatchMap must be sorted first.');
    }

    const self = this;

    return this.orderedResults.slice(0, n).map(result => {
      return new MatchResult(self, result.n, result.s);
    });
  }

  public bottom(n: number) {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    if (!this.orderedResults.length) {
      throw new Error('MatchMap must be sorted first.');
    }

    const self = this;

    return this.orderedResults
      .slice(this.orderedResults.length - n, n)
      .map(result => new MatchResult(self, result.n, result.s));
  }

  /* Can be optimized with binary splitting */

  public matchFrequencyData() {
    if (!this.initialized) {
      throw new Error('MatchMap must be initialized first.');
    }

    if (!this.orderedResults.length) {
      throw new Error('MatchMap must be sorted first.');
    }

    if (this.myMatchFrequencyData) {
      return this.myMatchFrequencyData;
    }

    const ordered = this.orderedResults;
    const matchFrequencyData = makeArray(this.query.size + 1);

    let maxMatch = this.query.size;
    let lastIndex = 0;
    let num;

    for (let i = 0, len = ordered.length; i < len; i++) {
      num = ordered[i].s;
      if (num < maxMatch) {
        matchFrequencyData[maxMatch] = i - lastIndex;
        lastIndex = i;
        maxMatch = num;
      }
      if (num === 0) {
        matchFrequencyData[0] = len - i;
        break;
      }
    }

    this.myMatchFrequencyData = matchFrequencyData;
    return this.myMatchFrequencyData;
  }

  public countMatches(int: number, aBitCount: Uint8Array) {
    int |= int >>> 1;
    int |= int >>> 2;
    int &= 0x11111111;
    int |= int >>> 3;
    int |= int >>> 6;
    return aBitCount[((int >>> 12) & 0xf0) | (int & 0xf)];
  }

  public execute(queryBuffer: any, searchSpaceBuffer: any) {
    let queryInts;
    let spaceInts;
    let queryIntsLength;
    let spaceIntsLength;
    let arrLen;
    let mapBuffer;
    let mapArray;
    let A;
    let B;
    let A1;
    let A2;
    let T;
    let cur;
    let pos;
    let i;
    let k;
    let adjustNeg;
    let adjustPos;
    let fnCountMatches;
    const bitCount = makeBitCount();

    queryInts = new Uint32Array(queryBuffer, 4);
    spaceInts = new Uint32Array(searchSpaceBuffer, 4);

    fnCountMatches = this.countMatches;

    queryIntsLength = queryInts.length | 0;
    spaceIntsLength = spaceInts.length | 0;

    arrLen = (queryIntsLength + spaceIntsLength) << 3;
    mapBuffer = new ArrayBuffer(4 * arrLen);
    mapArray = new Uint32Array(mapBuffer);

    for (k = 0 | 0; k < queryIntsLength; k++) {
      A = queryInts[k];
      cur = (queryIntsLength - k) << 3;

      for (i = 0 | 0; i < spaceIntsLength; i++) {
        // tslint:disable-next-line:no-unused-expression
        (T = A & spaceInts[i]) && (mapArray[(i << 3) + cur] += fnCountMatches(T, bitCount));
      }

      A1 = A >>> 4;
      A2 = A << 4;

      adjustNeg = cur - 1;
      adjustPos = cur + 1;

      while (A1 || A2) {
        for (i = 0 | 0; i < spaceIntsLength; i++) {
          B = spaceInts[i];
          pos = i << 3;
          // tslint:disable-next-line:no-unused-expression
          (T = A1 & B) && (mapArray[pos + adjustNeg] += fnCountMatches(T, bitCount));
          // tslint:disable-next-line:no-unused-expression
          (T = A2 & B) && (mapArray[pos + adjustPos] += fnCountMatches(T, bitCount));
        }

        A1 >>>= 4;
        A2 <<= 4;

        --adjustNeg;
        ++adjustPos;
      }
    }

    return mapBuffer;
  }
}

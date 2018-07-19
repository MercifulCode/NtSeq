import { MatchResult } from './MatchResult';
import { Seq } from './Seq';
export declare class MatchMap {
    readonly query: Seq;
    readonly searchSpace: Seq;
    myResults?: number[] | undefined;
    orderedResults: any[];
    myMatchFrequencyData: any[] | null;
    initialized: boolean;
    offset: number;
    positionAdjustment: number;
    debug: {
        prepareTime: number | null;
        searchTime: number | null;
        sortTime: number | null;
    };
    constructor(query: Seq, searchSpace: Seq, offset?: any);
    initialize(results?: number[]): this;
    sort(): this;
    __calculate_p_match(query: Seq, searchSpace: Seq): number;
    results(offset: number, count: number): number[];
    best(): MatchResult;
    top(n: number): MatchResult[];
    bottom(n: number): MatchResult[];
    matchFrequencyData(): any[];
    countMatches(int: number, aBitCount: Uint8Array): number;
    execute(queryBuffer: any, searchSpaceBuffer: any): ArrayBuffer;
}

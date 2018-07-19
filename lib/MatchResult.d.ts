import { MatchMap } from './MatchMap';
import { Seq } from './Seq';
export declare class MatchResult {
    readonly matchMap: MatchMap;
    readonly position: number;
    readonly score: number;
    align: any;
    s: number;
    constructor(matchMap: MatchMap, position: number, score: number);
    alignment(): any;
    alignmentMask(): Seq;
    alignmentCover(): Seq;
}

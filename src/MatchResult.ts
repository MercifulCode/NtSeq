import { MatchMap } from './MatchMap';
import { Seq } from './Seq';

// tslint:disable:no-bitwise
export class MatchResult {
  public align = null as any;
  public s = 0;
  constructor(readonly matchMap: MatchMap, readonly position: number, readonly score: number) {}

  public alignment() {
    if (!this.align) {
      const map = this.matchMap;
      if (this.position < 0) {
        this.align = new Seq()
          .read('-')
          .repeat(-this.position)
          .polymerize(map.searchSpace!.replicate(0, map.query!.length + this.position));
      } else if (this.position + map.query!.length > map.searchSpace!.length) {
        this.align = map
          .searchSpace!.replicate(this.position)
          .polymerize(new Seq().read('-').repeat(map.searchSpace!.length - this.position));
      } else {
        this.align = map.searchSpace!.replicate(this.position, map.query!.length);
      }
    }
    return this.align;
  }

  public alignmentMask() {
    return this.matchMap.query!.mask(this.alignment());
  }

  public alignmentCover() {
    return this.matchMap.query!.cover(this.alignment());
  }
}

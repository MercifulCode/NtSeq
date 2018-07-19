"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var Seq_1 = require("./Seq");
// tslint:disable:no-bitwise
var MatchResult = /** @class */ (function () {
    function MatchResult(matchMap, position, score) {
        this.matchMap = matchMap;
        this.position = position;
        this.score = score;
        this.align = null;
        this.s = 0;
    }
    MatchResult.prototype.alignment = function () {
        if (!this.align) {
            var map = this.matchMap;
            if (this.position < 0) {
                this.align = new Seq_1.Seq()
                    .read('-')
                    .repeat(-this.position)
                    .polymerize(map.searchSpace.replicate(0, map.query.length + this.position));
            }
            else if (this.position + map.query.length > map.searchSpace.length) {
                this.align = map
                    .searchSpace.replicate(this.position)
                    .polymerize(new Seq_1.Seq().read('-').repeat(map.searchSpace.length - this.position));
            }
            else {
                this.align = map.searchSpace.replicate(this.position, map.query.length);
            }
        }
        return this.align;
    };
    MatchResult.prototype.alignmentMask = function () {
        return this.matchMap.query.mask(this.alignment());
    };
    MatchResult.prototype.alignmentCover = function () {
        return this.matchMap.query.cover(this.alignment());
    };
    return MatchResult;
}());
exports.MatchResult = MatchResult;
//# sourceMappingURL=MatchResult.js.map
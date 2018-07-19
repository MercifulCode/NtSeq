"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var MatchResult_1 = require("./MatchResult");
var nt_1 = require("./nt");
// tslint:disable:no-bitwise
var MatchMap = /** @class */ (function () {
    function MatchMap(query, searchSpace, offset) {
        this.query = query;
        this.searchSpace = searchSpace;
        this.myResults = new Array();
        this.orderedResults = new Array();
        this.myMatchFrequencyData = null;
        this.initialized = false;
        this.offset = 0;
        this.positionAdjustment = 0;
        this.debug = {
            prepareTime: null,
            searchTime: null,
            sortTime: null,
        };
    }
    MatchMap.prototype.initialize = function (results) {
        this.orderedResults = [];
        this.myMatchFrequencyData = null;
        this.initialized = true;
        var t = new Date().valueOf();
        if (!results) {
            var dataArray = new Uint32Array(this.execute(this.query.buffer, this.searchSpace.buffer));
            this.debug.searchTime = -t + (t = new Date().valueOf());
            var queryLen = this.query.size;
            var searchLen = this.searchSpace.size;
            var resultsLen = ((8 - (queryLen % 8)) % 8) + 1;
            var totalLen = searchLen + queryLen - 1;
            results = [].slice.call(dataArray, resultsLen, resultsLen + totalLen);
        }
        this.myResults = results;
        this.debug.prepareTime = -t + (t = new Date().valueOf());
        return this;
    };
    MatchMap.prototype.sort = function () {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        var t = new Date().valueOf();
        if (this.debug.sortTime !== null) {
            return this;
        }
        var adjust = this.positionAdjustment;
        this.orderedResults = this.myResults
            .map(function (v, i) { return { n: i + adjust, s: v }; })
            .sort(function (a, b) { return b.s - a.s; });
        this.debug.sortTime = new Date().valueOf() - t;
        return this;
    };
    MatchMap.prototype.__calculate_p_match = function (query, searchSpace) {
        /*
        The approximate probability that two randomly chosen nucleotides
        from QUERY and SEARCH match each other
      */
        var queryContent = query.fractionalContentATGC();
        var searchSpaceContent = searchSpace.fractionalContentATGC();
        return (queryContent.A * searchSpaceContent.A +
            queryContent.T * searchSpaceContent.T +
            queryContent.G * searchSpaceContent.G +
            queryContent.C * searchSpaceContent.C);
    };
    MatchMap.prototype.results = function (offset, count) {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        if (offset === undefined) {
            return this.myResults.slice();
        }
        if (count === undefined) {
            return this.myResults.slice(offset | 0);
        }
        return this.myResults.slice(offset | 0, count | 0);
    };
    MatchMap.prototype.best = function () {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        if (!this.orderedResults.length) {
            throw new Error('MatchMap must be sorted first.');
        }
        var result = this.orderedResults[0];
        return new MatchResult_1.MatchResult(this, result.n, result.s);
    };
    MatchMap.prototype.top = function (n) {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        if (!this.orderedResults.length) {
            throw new Error('MatchMap must be sorted first.');
        }
        var self = this;
        return this.orderedResults.slice(0, n).map(function (result) {
            return new MatchResult_1.MatchResult(self, result.n, result.s);
        });
    };
    MatchMap.prototype.bottom = function (n) {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        if (!this.orderedResults.length) {
            throw new Error('MatchMap must be sorted first.');
        }
        var self = this;
        return this.orderedResults
            .slice(this.orderedResults.length - n, n)
            .map(function (result) { return new MatchResult_1.MatchResult(self, result.n, result.s); });
    };
    /* Can be optimized with binary splitting */
    MatchMap.prototype.matchFrequencyData = function () {
        if (!this.initialized) {
            throw new Error('MatchMap must be initialized first.');
        }
        if (!this.orderedResults.length) {
            throw new Error('MatchMap must be sorted first.');
        }
        if (this.myMatchFrequencyData) {
            return this.myMatchFrequencyData;
        }
        var ordered = this.orderedResults;
        var matchFrequencyData = nt_1.makeArray(this.query.size + 1);
        var maxMatch = this.query.size;
        var lastIndex = 0;
        var num;
        for (var i = 0, len = ordered.length; i < len; i++) {
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
    };
    MatchMap.prototype.countMatches = function (int, aBitCount) {
        int |= int >>> 1;
        int |= int >>> 2;
        int &= 0x11111111;
        int |= int >>> 3;
        int |= int >>> 6;
        return aBitCount[((int >>> 12) & 0xf0) | (int & 0xf)];
    };
    MatchMap.prototype.execute = function (queryBuffer, searchSpaceBuffer) {
        var queryInts;
        var spaceInts;
        var queryIntsLength;
        var spaceIntsLength;
        var arrLen;
        var mapBuffer;
        var mapArray;
        var A;
        var B;
        var A1;
        var A2;
        var T;
        var cur;
        var pos;
        var i;
        var k;
        var adjustNeg;
        var adjustPos;
        var fnCountMatches;
        var bitCount = nt_1.makeBitCount();
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
                    (T = A1 & B) && (mapArray[pos + adjustNeg] += fnCountMatches(T, bitCount));
                    (T = A2 & B) && (mapArray[pos + adjustPos] += fnCountMatches(T, bitCount));
                }
                A1 >>>= 4;
                A2 <<= 4;
                --adjustNeg;
                ++adjustPos;
            }
        }
        return mapBuffer;
    };
    return MatchMap;
}());
exports.MatchMap = MatchMap;
//# sourceMappingURL=MatchMap.js.map
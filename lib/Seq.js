"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
// tslint:disable:no-bitwise
var MatchMap_1 = require("./MatchMap");
var nt = require("./nt");
var fs = require("fs");
var Seq = /** @class */ (function () {
    function Seq(type) {
        if (type === void 0) { type = 'DNA'; }
        this.type = type;
        this.buffer = new ArrayBuffer(4);
        this.endPadding = 0;
        this.length = 0;
        this.myComplement = null;
        this.myContent = null;
        this.myContentATGC = null;
        this.myFractionalContent = null;
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
    Object.defineProperty(Seq.prototype, "isRNA", {
        get: function () {
            return this.type === 'RNA';
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Seq.prototype, "size", {
        get: function () {
            return this.length;
        },
        enumerable: true,
        configurable: true
    });
    Seq.prototype.read = function (strData) {
        var ntToByte = nt.nucleotidesToByte;
        var nucleotideString = strData
            .toUpperCase()
            .replace(/\s/g, '')
            .replace('U', 'T')
            .replace(/[^ATGCBVHDMKRYSWN\-]/g, '-');
        var length = nucleotideString.length | 0;
        var max = length >>> 1;
        var odd = length & 1;
        var endPadding = (4 - ((max + odd) % 4)) % 4;
        this.endPadding = endPadding;
        var buffer = new ArrayBuffer(max + odd + endPadding + 4);
        var dataArray = new Int8Array(buffer, 4);
        var n;
        var i = 0;
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
    };
    Seq.prototype.readBuffer = function (buffer) {
        this.buffer = buffer;
        var length = new Uint32Array(buffer, 0, 1)[0];
        var max = length >>> 1;
        var odd = length & 1;
        var endPadding = (4 - ((max + odd) % 4)) % 4;
        this.endPadding = endPadding;
        this.length = length;
        this.myComplement = null;
        this.myContent = null;
        this.myFractionalContent = null;
        this.myContentATGC = null;
        this.myFractionalContent = null;
        return this;
    };
    Seq.prototype.readFASTA = function (strFASTA) {
        var data = strFASTA.split(/\n\r?/gi);
        while (data.length && data[0][0] === '>') {
            data.shift();
        }
        return this.read(data.join(''));
    };
    Seq.prototype.byteComplement = function () {
        var bComp = nt.NT.byteComplement;
        var fwdBuffer = this.buffer;
        var len = fwdBuffer.byteLength;
        var n;
        var i;
        var copyBuffer;
        var fromArray;
        var copyArray;
        var isOdd = this.length & 1;
        if (isOdd) {
            copyBuffer = new ArrayBuffer(len);
            copyArray = new Uint32Array(copyBuffer, 4);
            fromArray = new Uint32Array(fwdBuffer, 4);
            n = (len - 4) >>> 2;
            while (n--) {
                copyArray[n] = (fromArray[n] << 4) | (fromArray[n - 1] >>> 28);
            }
        }
        else {
            copyBuffer = fwdBuffer;
        }
        var fwdArray = new Uint8Array(copyBuffer, 4);
        var buffer = new ArrayBuffer(len);
        var dataArray = new Uint8Array(buffer, 4);
        n = len - 4 - this.endPadding;
        i = 0;
        while (n--) {
            dataArray[i++] = bComp[fwdArray[n]];
        }
        new Uint32Array(buffer, 0, 1)[0] = this.length;
        return buffer;
    };
    Seq.prototype.sequence = function () {
        var byteToNt = nt.NT.byteToNucleotides;
        var buffer = this.buffer;
        if (buffer.byteLength < 4) {
            return '';
        }
        var dataArray = new Uint8Array(buffer, 4);
        var len = buffer.byteLength - 4 - this.endPadding;
        var nts = nt.makeArray(len);
        var i = 0;
        do {
            nts[i] = byteToNt[dataArray[i]];
            ++i;
        } while (i < len);
        var returnString;
        i = nts.length - 1;
        if (this.length & 1) {
            nts[i] = nts[i][0];
        }
        returnString = nts.join('');
        if (this.isRNA) {
            returnString = returnString.replace(/T/gi, 'U');
        }
        return returnString;
    };
    Object.defineProperty(Seq.prototype, "complement", {
        get: function () {
            if (!this.myComplement) {
                this.myComplement = this.byteComplement();
            }
            var complement = new Seq(this.type).readBuffer(this.myComplement.slice(0));
            complement.myComplement = this.buffer.slice(0);
            return complement;
        },
        enumerable: true,
        configurable: true
    });
    Seq.prototype.equivalent = function (seq) {
        if (!(seq instanceof Seq)) {
            throw new Error('Can only check for equivalence between sequences');
        }
        if (this.type !== seq.type) {
            return false;
        }
        var checkInts = new Uint32Array(this.buffer);
        var compareInts = new Uint32Array(seq.buffer);
        for (var i = 0, len = checkInts.length; i < len; i++) {
            if (checkInts[i] !== compareInts[i]) {
                return false;
            }
        }
        return true;
    };
    Seq.prototype.replicate = function (start, length) {
        if (start === void 0) { start = 0; }
        start |= 0;
        if (start < 0) {
            start = this.length + start;
        }
        if (length === undefined) {
            if (start === 0) {
                return this.clone();
            }
            length = this.length - start;
        }
        else {
            length |= 0;
            length = Math.min(length, this.length - start);
        }
        length = Math.min(length, this.length - start);
        if (length <= 0) {
            return this.nullSeq();
        }
        return this.slice(start, length);
    };
    Seq.prototype.polymerize = function (seq) {
        var seqLen = seq.length;
        if (!(seq instanceof Seq)) {
            throw new Error('.polymerize requires valid sequence');
        }
        if (!this.length) {
            return seq.clone();
        }
        var length = this.length + seqLen;
        var max = length >>> 1;
        var odd = length & 1;
        var endPadding = (4 - ((max + odd) % 4)) % 4;
        var newBuffer = new ArrayBuffer(max + odd + endPadding + 4);
        var newArray = new Uint32Array(newBuffer, 4);
        var copyBuffer = this.buffer;
        var copyArray = new Uint32Array(copyBuffer, 4);
        var seqBuffer = seq.buffer;
        var seqArray = new Uint32Array(seqBuffer, 4);
        var copyPos = 0;
        var shift = (this.length % 8) * 4;
        var shiftSeq = 32 - shift;
        for (var len = copyArray.length; copyPos < len; copyPos++) {
            newArray[copyPos] = copyArray[copyPos];
        }
        if (shift) {
            newArray[--copyPos] |= seqArray[0] << shift;
            for (var i = 0, len = seqArray.length; i < len; i++) {
                newArray[++copyPos] = (seqArray[i] >>> shiftSeq) | (seqArray[i + 1] << shift);
            }
        }
        else {
            for (var i = 0, len = seqArray.length; i < len; i++) {
                newArray[copyPos++] = seqArray[i];
            }
        }
        new Uint32Array(newBuffer, 0, 1)[0] = length;
        return new Seq(this.type).readBuffer(newBuffer);
    };
    Seq.prototype.insertion = function (seq, offset) {
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
    };
    Seq.prototype.deletion = function (offset, count) {
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
    };
    Seq.prototype.repeat = function (count) {
        count |= 0;
        var copy = this.replicate();
        var base = new Seq(this.type);
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
    };
    Seq.prototype.mask = function (seq) {
        if (!(seq instanceof Seq)) {
            throw new Error('Can only mask with valid sequence');
        }
        var newBuffer = this.buffer.slice(0);
        var newArray = new Uint32Array(newBuffer, 4);
        var compareArray = new Uint32Array(seq.buffer, 4);
        for (var i = 0, len = newArray.length; i < len; i++) {
            newArray[i] &= compareArray[i];
        }
        return new Seq(this.type).readBuffer(newBuffer);
    };
    Seq.prototype.cover = function (seq) {
        if (!(seq instanceof Seq)) {
            throw new Error('Can only cover with valid sequence');
        }
        var newBuffer = this.buffer.slice(0);
        var newArray = new Uint32Array(newBuffer, 4);
        var compareArray = new Uint32Array(seq.buffer, 4);
        for (var i = 0, len = newArray.length; i < len; i++) {
            newArray[i] |= compareArray[i];
        }
        return new Seq(this.type).readBuffer(newBuffer);
    };
    Seq.prototype.nullSeq = function () {
        return new Seq(this.type).readBuffer(new ArrayBuffer(4));
    };
    Seq.prototype.clone = function () {
        return new Seq(this.type).readBuffer(this.buffer.slice(0));
    };
    Seq.prototype.slice = function (start, length) {
        var max = length >>> 1;
        var odd = length & 1;
        var endPadding = (4 - ((max + odd) % 4)) % 4;
        var newBuffer = new ArrayBuffer(max + odd + endPadding + 4);
        var newArray = new Uint32Array(newBuffer, 4);
        var subBuffer = this.buffer.slice(4 + (start >>> 1), 4 + (start >>> 1) + newBuffer.byteLength);
        var subInt32Length = subBuffer.byteLength >>> 2;
        var subArray = new Uint32Array(subBuffer, 0, subInt32Length);
        if (start & 1) {
            var subArrayIndex = 0;
            do {
                newArray[subArrayIndex] = (subArray[subArrayIndex] >>> 4) | (subArray[subArrayIndex + 1] << 28);
                ++subArrayIndex;
            } while (subArrayIndex < subArray.length);
            var remainder = subBuffer.byteLength - subArray.byteLength;
            if (remainder) {
                var remainderArray = new Uint8Array(newBuffer, 4 + (subArrayIndex << 2));
                var subRemainderArray = new Uint8Array(subBuffer, subArrayIndex << 2);
                if (newArray.length > 0) {
                    newArray[subArrayIndex - 1] |= subRemainderArray[0] << 28;
                }
                for (var i = 0, len = subRemainderArray.length; i < len; i++) {
                    remainderArray[i] = (subRemainderArray[i] >>> 4) | (subRemainderArray[i + 1] << 4);
                }
            }
        }
        else {
            var subArrayIndex = 0;
            do {
                newArray[subArrayIndex] = subArray[subArrayIndex];
                ++subArrayIndex;
            } while (subArrayIndex < subArray.length);
            var remainder = subArray.byteLength - subBuffer.byteLength;
            if (remainder) {
                var remainderArray = new Uint8Array(newBuffer, 4 + (subArrayIndex << 2));
                var subRemainderArray = new Uint8Array(subBuffer, subArrayIndex << 2);
                for (var i = 0, len = subRemainderArray.length; i < len; i++) {
                    remainderArray[i] = subRemainderArray[i];
                }
            }
        }
        var clearShift = (endPadding * 2 + odd) * 4;
        var clearOut = new Uint32Array(newBuffer, newBuffer.byteLength - 4);
        clearOut[0] = (clearOut[0] << clearShift) >>> clearShift;
        new Uint32Array(newBuffer, 0, 1)[0] = length;
        return new Seq(this.type).readBuffer(newBuffer);
    };
    Object.defineProperty(Seq.prototype, "content", {
        get: function () {
            if (!this.myContent) {
                var ntContentByte = nt.makeArray(256);
                var buffer = this.buffer;
                var dataArray = new Uint8Array(buffer);
                for (var i = 4; i < buffer.byteLength - this.endPadding; i++) {
                    ntContentByte[dataArray[i]]++;
                }
                var binToNt = nt.binToNucleotide;
                var ntList = nt.nucleotideList;
                var ntContent = Object.create(null);
                for (var i = 0, len = ntList.length; i < len; i++) {
                    ntContent[ntList[i]] = 0;
                }
                for (var i = 0, len = ntContentByte.length; i < len; i++) {
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
            var returnContent = Object.create(null);
            var keys = Object.keys(this.myContent);
            for (var i = 0, len = keys.length; i < len; i++) {
                returnContent[keys[i]] = this.myContent[keys[i]];
            }
            return returnContent;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Seq.prototype, "fractionalContent", {
        get: function () {
            if (!this.myFractionalContent) {
                var content = this.content;
                var nts = Object.keys(content);
                for (var i = 0, len = nts.length; i < len; i++) {
                    content[nts[i]] = content[nts[i]] / this.length;
                }
                this.myFractionalContent = content;
            }
            var returnContent = Object.create(null);
            var keys = Object.keys(this.myFractionalContent);
            for (var i = 0, len = keys.length; i < len; i++) {
                returnContent[keys[i]] = this.myFractionalContent[keys[i]];
            }
            return returnContent;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Seq.prototype, "contentATGC", {
        get: function () {
            if (!this.myContentATGC) {
                var ntToBin = nt.nucleotideToBin;
                var content = this.content;
                var nts = Object.keys(content);
                var contentATGC = Object.create(null);
                contentATGC.A = 0;
                contentATGC.T = 0;
                contentATGC.G = 0;
                contentATGC.C = 0;
                var bits = 0;
                var nucleotides = void 0;
                var ntBin = void 0;
                var n = void 0;
                var curContent = void 0;
                for (var i = 0, len = nts.length; i < len; i++) {
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
            var returnContent = Object.create(null);
            var keys = Object.keys(this.myContentATGC);
            for (var i = 0, len = keys.length; i < len; i++) {
                returnContent[keys[i]] = this.myContentATGC[keys[i]];
            }
            return returnContent;
        },
        enumerable: true,
        configurable: true
    });
    Object.defineProperty(Seq.prototype, "fractionalContentATGC", {
        get: function () {
            if (!this.myFractionalContentATGC) {
                var content = this.contentATGC;
                var nts = Object.keys(content);
                for (var i = 0, len = nts.length; i < len; i++) {
                    content[nts[i]] = content[nts[i]] / this.length;
                }
                this.myFractionalContentATGC = content;
            }
            var returnContent = Object.create(null);
            var keys = Object.keys(this.myFractionalContentATGC);
            for (var i = 0, len = keys.length; i < len; i++) {
                returnContent[keys[i]] = this.myFractionalContentATGC[keys[i]];
            }
            return returnContent;
        },
        enumerable: true,
        configurable: true
    });
    Seq.prototype.translate = function (ntOffset, ntCount) {
        if (ntOffset === void 0) { ntOffset = 0; }
        var binToAA = nt.__12BitToAminoAcid;
        ntOffset |= 0;
        if (ntCount === undefined) {
            ntCount = this.length - ntOffset;
        }
        ntCount |= 0;
        ntCount -= ntCount % 3;
        var offset = (ntOffset >>> 1) + 4;
        var max = offset + (ntCount >>> 1) + (ntCount & 1);
        var dataArray = new Uint8Array(this.buffer);
        var aminoAcids = nt.makeArray(ntCount / 3);
        /**/
        var aa = 0;
        var lastByte;
        var byte1;
        var byte2;
        var byte3;
        if ((ntOffset & 1) === 0) {
            for (var i = offset; i < max; i += 3) {
                byte1 = dataArray[i];
                byte2 = dataArray[i + 1];
                byte3 = dataArray[i + 2];
                aminoAcids[aa++] = binToAA[byte1 | ((byte2 & 0xf) << 8)];
                aminoAcids[aa++] = binToAA[(byte3 << 4) | (byte2 >>> 4)];
            }
        }
        else {
            lastByte = dataArray[offset];
            for (var i = offset + 1; i < max; i += 3) {
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
    };
    Seq.prototype.translateFrame = function (frame, AAoffset, AAcount) {
        if (frame === void 0) { frame = 0; }
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
    };
    Seq.prototype.mapSequence = function (seq, offset) {
        if (!(seq instanceof Seq)) {
            throw new Error('.mapSequence requires valid Seq');
        }
        return new MatchMap_1.MatchMap(seq, this, offset);
    };
    Seq.prototype.loadFile = function (path, ext) {
        if (ext === 'fasta') {
            return this.loadFASTA(path);
        }
        else {
            return this.load4bnt(path);
        }
    };
    ;
    Seq.prototype.loadFASTA = function (path) {
        return this.readFASTA(fs.readFileSync(path).toString());
    };
    ;
    Seq.prototype.load4bnt = function (path) {
        var nodeBuffer = fs.readFileSync(path);
        var buffer = new ArrayBuffer(nodeBuffer.length);
        var view = new Uint8Array(buffer);
        for (var i = 0; i < nodeBuffer.length; ++i) {
            view[i] = nodeBuffer[i];
        }
        return this.readBuffer(buffer);
    };
    return Seq;
}());
exports.Seq = Seq;
//# sourceMappingURL=Seq.js.map
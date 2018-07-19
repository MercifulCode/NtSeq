import { MatchMap } from '../MatchMap';
import { Seq } from '../Seq';

const naiveSearch = (query: string, space: string) => {
  var qLen = query.length;
  var sLen = space.length;
  var map = new Uint32Array(qLen + sLen);

  var curChar;
  var offset;

  for(var j = 0; j < qLen; j++) {
    curChar = query[j];
    offset = qLen - j;
    for(var i = 0; i < sLen; i++) {
      (curChar === space[i]) && (++map[offset + i]);
    }
  }
  return map;
}

const inputs = [

  {name: '1,000,000, 0%', data: [26, 2501, 'TTTT', 'AAAA']},
  {name: '10,000,000, 0%', data: [26, 25001, 'TTTT', 'AAAA']},
  {name: '100,000,000, 0%', data: [26, 250001, 'TTTT', 'AAAA']},
  //{name: '1,000,000,000, 0%', data: [26, 2500001, 'TTTT', 'AAAA']},

  {name: '1,000,000, 25%', data: [26, 2501, 'ATGC', 'AAAA']},
  {name: '10,000,000, 25%', data: [26, 25001, 'ATGC', 'AAAA']},
  {name: '100,000,000, 25%', data: [26, 250001, 'ATGC', 'AAAA']},
  //{name: '1,000,000,000, 25%', data: [26, 2500001, 'ATGC', 'AAAA']},

  {name: '1,000,000, 50%', data: [26, 2501, 'ATAT', 'AAAA']},
  {name: '10,000,000, 50%', data: [26, 25001, 'ATAT', 'AAAA']},
  {name: '100,000,000, 50%', data: [26, 250001, 'ATAT', 'AAAA']},
  //{name: '1,000,000,000, 50%', data: [26, 2500001, 'ATAT', 'AAAA']},

  {name: '1,000,000, 100%', data: [26, 2501, 'AAAA', 'AAAA']},
  {name: '10,000,000, 100%', data: [26, 25001, 'AAAA', 'AAAA']},
  {name: '100,000,000, 100%', data: [26, 250001, 'AAAA', 'AAAA']},
  //{name: '1,000,000,000, 100%', data: [26, 2500001, 'AAAA', 'AAAA']},

];

describe('Benchmarking', () => {
  it('Should benchmark.', () => {
    var results: {[key: string]: any} = new Array();

    for(var i = 0, len = inputs.length; i < len; i++) {

      var qString = Array(inputs[i].data[0]).join(inputs[i].data[2].toString());
      var sString = Array(inputs[i].data[1]).join(inputs[i].data[3].toString());
      var a = new Seq().read(qString);
      var b = new Seq().read(sString);

      console.log('Benchmark ' + (i + 1) + ' of ' + len + '...')

      var naiveTime = (new Date).valueOf();
      naiveSearch(qString, sString);
      naiveTime = (new Date).valueOf() - naiveTime;
      console.log('\t ... Naive search complete!')

      var map = new MatchMap(a, b);
      console.log('\t ... Nt.MatchMap complete!')

      results.push({
        name: inputs[i].name,
        search: map.debug.searchTime,
        searchScore: ((map.debug.searchTime! * 1000000) / Number.parseInt((a.size * b.size).toFixed(2))),
        naive: naiveTime,
        naiveScore: ((naiveTime * 1000000) / Number.parseInt((a.size * b.size).toFixed(2))),
        prepare: map.debug.prepareTime,
        prepareScore: ((map.debug.prepareTime! * 1000000) / Number.parseInt((a.size * b.size).toFixed(2))),
        sort: map.debug.sortTime,
        sortScore: ((map.debug.sortTime! * 1000000) / Number.parseInt((a.size * b.size).toFixed(2))),
      });

    }

    console.log('\n\n');

    var nameCellSize = 20;
    var cellSize = 12;

    var items = [
      ['naive', 'ms'],
      ['search', 'ms'],
      ['prepare', 'ms'],
      ['sort', 'ms'],
      ['naiveScore', 'ns/nt'],
      ['searchScore', 'ns/nt'],
      ['prepareScore', 'ns/nt'],
      ['sortScore', 'ns/nt'],
    ];

    var strArray = [];
    strArray.push(Array(nameCellSize - 'Benchmark'.length + 1).join(' ') + 'Benchmark');

    for(var j = 0, jlen = items.length; j < jlen; j++) {
      val = items[j][0];
      if(val.length > cellSize) {
        val = val.substring(0, cellSize - 3) + '...';
      }
      strArray.push(Array(cellSize - val.length + 1).join(' ') + val);
    }

    console.log(strArray.join(' | '));
    console.log(Array(strArray.join(' | ').length + 1).join('-'));

    for (let i = 0, len = results.length; i < len; i++) {

      var val;
      var result = results[i];

      var strArray = [];
      strArray.push(Array(nameCellSize - result.name.length + 1).join(' ') + result.name);

      for(let j = 0, jlen = items.length; j < jlen; j++) {
        if (result[items[j][0]]) {
          val = result[items[j][0]].toString() + items[j][1];
          strArray.push(Array(cellSize - val.length + 1).join(' ') + val);
        }
      }

      console.log(strArray.join(' | '));

    }

  })
})




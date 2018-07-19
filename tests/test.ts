var Nt = require('../lib/nt.js');
var test = require('./test_object.js');


test.add(function() {

  var a = new Nt.Seq().loadFASTA(__dirname + '/data/sequence.fasta');
  var b = new Nt.Seq().read('TCTTATTTGTGCTGTTTATT');

  // Sequence target in Bacteriophage T4, gene 46
  // [33302..34984]

  var matchMap = a.mapSequence(b).initialize().sort();

  var testData = [
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

  var isCorrect = true;
  var matchFrequencyData = matchMap.matchFrequencyData();

  for (var i = 0; i < matchFrequencyData.length; i++) {
    if (matchFrequencyData[i] !== testData[i]) {
      isCorrect = false;
      break;
    }
  }

  return isCorrect &&
    matchMap.best().position === 34767 &&
    matchMap.best().score === 20;

}, 'mapSequence successfully found correct match scores, loaded from .fasta');

test.add(function() {

  var a = new Nt.Seq().load4bnt(__dirname + '/data/sequence.4bnt');
  var b = new Nt.Seq().read('TCTTATTTGTGCTGTTTATT');

  // Sequence target in Bacteriophage T4, gene 46
  // [33302..34984]

  var matchMap = a.mapSequence(b).initialize().sort();

  var testData = [
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

  var isCorrect = true;
  var matchFrequencyData = matchMap.matchFrequencyData();

  for (var i = 0; i < matchFrequencyData.length; i++) {
    if (matchFrequencyData[i] !== testData[i]) {
      isCorrect = false;
      break;
    }
  }

  return isCorrect &&
    matchMap.best().position === 34767 &&
    matchMap.best().score === 20;

}, 'mapSequence successfully found correct match scores, loaded from .4bnt');

test.run();

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import logging
import numpy as np

from pbcore.io import BasH5Reader, BaxH5Reader
from pbcore.io.BarcodeH5Reader import LabeledZmw
from pbbarcode.SWaligner import SWaligner
from pbbarcode.utils import makeBarcodeLabel, Bunch, reverseComplement

def makeFromRangeFunc(barcodeLength, insertSidePad, adapterSidePad, useOldWorkflow):
    """In order to score an Adapter for possible barcodes, we need a function that
    returns the ranges of sequence immediately to the 5' and 3' end of a given
    adapter, which differ slightly between workflows."""

    if useOldWorkflow:
        # The old fromRange function reports ranges in their default orientation
        def fromRangeFunc(zmw, rStart, rEnd):
            try:
                qSeqLeft = zmw.read(rStart - (barcodeLength + insertSidePad),
                                    rStart + adapterSidePad).basecalls()
            except IndexError:
                qSeqLeft = None
            try:
                qSeqRight = zmw.read(rEnd - adapterSidePad,
                                     rEnd + barcodeLength + insertSidePad).basecalls()
            except IndexError:
                qSeqRight = None
            return (qSeqLeft, qSeqRight)
    else:
        # The new fromRange function reports ranges oriented away from the Adapter
        def fromRangeFunc(zmw, rStart, rEnd):
            try:
                qSeqLeftRaw = zmw.read(rStart - (barcodeLength + insertSidePad),
                                       rStart + adapterSidePad).basecalls()
                qSeqLeft = reverseComplement( qSeqLeftRaw )
            except IndexError:
                qSeqLeft = None
            try:
                qSeqRight = zmw.read(rEnd - adapterSidePad,
                                     rEnd + barcodeLength + insertSidePad).basecalls()
            except IndexError:
                qSeqRight = None
            return (qSeqLeft, qSeqRight)

    # Return the selected fromRange function
    return fromRangeFunc

def makeScoreFlankingFunc(forwardScorer, reverseScorer, pairedScorer,
                          scoreMode, numSeqs, oldWorkflow):
    """Once the flanking regions around an adapter have been extracted
    they need to be scored against the appropriate set of barcode sequences,
    the specifics of which vary by scoreMode and workflow.
    """

    def scoreFlankingOld(holeNum, adapters, scoredFirst):
        adapterScores = [[]]*len(adapters)
        barcodeScores = np.zeros(numSeqs)
        for i, adapter in enumerate(adapters):
            fscores  = forwardScorer(adapter[0])
            rscores  = reverseScorer(adapter[0])
            ffscores = forwardScorer(adapter[1])
            rrscores = reverseScorer(adapter[1])

            # Average the two flanking scores for the adapter score
            if adapter[0] and adapter[1]:
                adapterScores[i] = np.maximum((fscores + rrscores)/2.0,
                                             (rscores + ffscores)/2.0)
            # Single flanking scores are taken as-is
            elif adapter[0] or adapter[1]:
                adapterScores[i] = np.maximum((fscores + rrscores),
                                             (rscores + ffscores))
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else barcodeScores

        return Bunch(holeNum=holeNum, numAdapters=len(adapters),
                     barcodeScores=barcodeScores, adapterScores=adapterScores,
                     scoredFirst=scoredFirst)

    def scoreFlankingPaired(holeNum, adapters, scoredFirst):
        adapterScores = [[]]*len(adapters)
        barcodeScores = np.zeros(numSeqs)
        for i, adapter in enumerate(adapters):
            fscores = pairedScorer(adapter[0])
            rscores = pairedScorer(adapter[1])

            # Average the two flanking scores for the adapter score
            if adapter[0] and adapter[1]:
                adapterScores[i] = (fscores + rscores)/2.0
            # Single flanking scores are taken as-is
            elif adapter[0] or adapter[1]:
                adapterScores[i] = (fscores + rscores)
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else barcodeScores

        return Bunch(holeNum=holeNum, numAdapters=len(adapters),
                     barcodeScores=barcodeScores, adapterScores=adapterScores,
                     scoredFirst=scoredFirst)

    def scoreFlankingSymmetric(holeNum, adapters, scoredFirst):
        adapterScores = [[]]*len(adapters)
        barcodeScores = np.zeros(numSeqs)
        for i, adapter in enumerate(adapters):
            fscores = forwardScorer(adapter[0])
            rscores = forwardScorer(adapter[1])

            # Average the two flanking scores for the adapter score
            if adapter[0] and adapter[1]:
                adapterScores[i] = (fscores + rscores)/2.0
            # Single flanking scores are taken as-is
            elif adapter[0] or adapter[1]:
                adapterScores[i] = (fscores + rscores)
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else barcodeScores

        return Bunch(holeNum=holeNum, numAdapters=len(adapters),
                     barcodeScores=barcodeScores, adapterScores=adapterScores,
                     scoredFirst=scoredFirst)

    # Return the workflow/scoreMode appropriate scoring function
    if oldWorkflow:
        return scoreFlankingOld
    elif scoreMode == 'paired' and not oldWorkflow:
        return scoreFlankingPaired
    elif scoreMode == 'symmetric' and not oldWorkflow:
        return scoreFlankingSymmetric

# The following two functs create labeledZmws from a scoreBunch object
def makeSymmetricZmw(scoreBunch):
    """Convert a dictionary-like object with barcode scoring information
    into a LabeledZmw for a symmetrically barcoded read"""
    rankedBarcodes = np.argsort(-1 * scoreBunch.barcodeScores)
    bestIdx = rankedBarcodes[0]
    secondBestIdx = rankedBarcodes[1]
    return LabeledZmw(scoreBunch.holeNum,
                      scoreBunch.numAdapters,
                      bestIdx,
                      scoreBunch.barcodeScores[bestIdx],
                      secondBestIdx,
                      scoreBunch.barcodeScores[secondBestIdx],
                      scoreBunch.adapterScores)

def makePairedZmw(scoreBunch):
    """Convert a dictionary-like object with barcode scoring information
    into a LabeledZmw for a symmetrically barcoded read"""
    numSeqs = len(scoreBunch.barcodeScores)
    if scoreBunch.numAdapters == 1:
        # If we have one adapter, pick the best barcode from each pair as the
        #    score for that pair
        rawPairScores = [max(scoreBunch.barcodeScores[i], scoreBunch.barcodeScores[i+1]) \
                         for i in xrange(0, numSeqs, 2)]
        pairScores = np.array(rawPairScores)
        barcodeRanks = np.argsort(-pairScores)
        pairScores = pairScores[barcodeRanks]
    else:
        # If we have more than one adapter, score the two possible orderings we
        #    expec (F--R--F... or R--F--R...) then take the best score from
        #    those two possibilities as the score for that pair.
        # NOTE: A missed adapter will confuse this computation.
        scores  = scoreBunch.adapterScores
        results = np.zeros(numSeqs/2)
        for i in xrange(0, numSeqs, 2):
            orientations = [0,0]
            for j in xrange(0, len(scores)):
                orientations[j % 2] += scores[j][i]
                orientations[1 - j % 2] += scores[j][i + 1]
            results[i/2] = max(orientations)
        barcodeRanks = np.argsort(-results)
        pairScores = results[barcodeRanks]

    bestIdx = barcodeRanks[0]
    bestScore = pairScores[0]
    secondBestIdx = barcodeRanks[1]
    secondBestScore = pairScores[1]
    return LabeledZmw(scoreBunch.holeNum,
                      scoreBunch.numAdapters,
                      bestIdx,
                      bestScore,
                      secondBestIdx,
                      secondBestScore,
                      scoreBunch.adapterScores)


class BarcodeScorer(object):
    """A BarcodeScorer object scores ZMWs and produces summaries
    of the scores. Various parameters control the behavior of the
    object, specifically the padding allows the user to add a
    little extra on each side of the adapter find for safety. The
    most relevant parameter is the scoreMode which dictates how
    the barcodes are scored, either paired or symmetric."""
    def __init__(self, basH5, barcodeFasta,
                 adapterSidePad = 0,
                 insertSidePad = 4,
                 scoreMode = 'symmetric',
                 maxHits = 10,
                 scoreFirst = False,
                 startTimeCutoff = 1,
                 useOldWorkflow = False):

        self.basH5           = basH5
        self.barcodeFasta    = list(barcodeFasta)
        self.numSeqs         = len(self.barcodeFasta)
        self.barcodeNames    = np.array([x.name for x in self.barcodeFasta])
        self.aligner         = SWaligner(useOldWorkflow)
        self.useOldWorkflow  = useOldWorkflow
        self.adapterSidePad  = adapterSidePad
        self.insertSidePad   = insertSidePad
        self.maxHits         = maxHits
        self.scoreFirst      = scoreFirst
        self.startTimeCutoff = startTimeCutoff

        if scoreMode not in ['symmetric', 'paired']:
            raise Exception("scoreMode must either be symmetric or paired")
        self.scoreMode = scoreMode

        barcodeLengths = np.unique(map(lambda x : len(x.sequence),
                                       self.barcodeFasta))
        if len(barcodeLengths) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(barcodeLengths)

        # Original barcode sequences and scorers for the Old Workflow
        self.barcodeSeqs = [(barcode.sequence.upper(),
                             reverseComplement(barcode.sequence.upper()))
                            for barcode in self.barcodeFasta]
        forwardScorer = self.aligner.makeScorer([x[0] for x in self.barcodeSeqs])
        reverseScorer = self.aligner.makeScorer([x[1] for x in self.barcodeSeqs])

        # Forward-oriented barcode sequence pairs for the New Workflow
        self.orientedSeqs  = [bc.sequence.upper() if (i%2) == 0 else
                              reverseComplement(bc.sequence.upper())
                              for i, bc in enumerate(self.barcodeFasta)]
        pairedScorer  = self.aligner.makeScorer( self.orientedSeqs )

        # Given the scoreMode, create all of the possible barcode labels
        if self.scoreMode == 'paired':
            self.barcodeLabels = np.array([makeBarcodeLabel(self.barcodeFasta[i].name,
                                           self.barcodeFasta[i+1].name)
                                           for i in xrange(0, len(self.barcodeSeqs), 2)])
        else:
            self.barcodeLabels = np.array([makeBarcodeLabel(x.name, x.name)
                                           for x in self.barcodeFasta])

        # Make a "fromRange" function for finding adapter-flanking regions
        self.fromRange = makeFromRangeFunc(self.barcodeLength,
                                           self.insertSidePad,
                                           self.adapterSidePad,
                                           self.useOldWorkflow)

        # Make a "scoreFlankingRegions" function for the results of "fromRange"
        self.scoreFlankingRegions = makeScoreFlankingFunc(forwardScorer,
                                                          reverseScorer,
                                                          pairedScorer,
                                                          self.scoreMode,
                                                          self.numSeqs,
                                                          self.useOldWorkflow)

        # Select the score-mode appropriate function for formatting scoring
        #    results into LabeledZmw objects
        if self.scoreMode == 'paired':
            self.makeLabeledZmw = makePairedZmw
        else:
            self.makeLabeledZmw = makeSymmetricZmw

        # If initialization made it this far, log the settings used
        logging.debug(("Constructed BarcodeScorer with scoreMode: %s," + \
                "adapterSidePad: %d, insertSidePad: %d, scoreFirst: %r, and oldWorkflow: %s") \
                % (scoreMode, adapterSidePad, insertSidePad, scoreFirst, useOldWorkflow))

    @property
    def movieName(self):
        return self.basH5.movieName

    def _flankingSeqs(self, zmw):
        """Extract the flanking sequences for the first 'maxHits' adapters from a
        ZMW.  If that number is 0 and scoreFirst is true, try to extract a barcode
        from the 5' tip of the read instead.
        """
        # Extract the first X adapters
        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]

        # Extract the (left, right) sequence pairs around each adapter
        seqs = [self.fromRange(zmw, start, end) for (start, end) in adapterRegions]

        # We only score the first barcode if we don't find any adapters
        # *and* the start time is less than the threshold.
        scoredFirst = False
        if self.scoreFirst and not len(seqs):
            s = zmw.zmwMetric('HQRegionStartTime')
            e = zmw.zmwMetric('HQRegionEndTime')
            # s<e => has HQ.
            if s < e and s <= self.startTimeCutoff:
                l = self.barcodeLength + self.insertSidePad
                l = l if zmw.hqRegion[1] > l else zmw.hqRegion[1]
                try:
                    bc = zmw.read(0, l).basecalls()
                    if len(bc) >= self.barcodeLength:
                        seqs.insert(0, (bc, None))
                        scoredFirst = True
                except IndexError:
                    pass

        return (seqs, scoredFirst)

    def scoreZmw(self, zmw):
        flankingRegions, scoredFirst = self._flankingSeqs(zmw)
        return self.scoreFlankingRegions(zmw.holeNumber, flankingRegions, scoredFirst)

    def labelZmws(self, holeNumbers):
        """Return a list of LabeledZmws for input holeNumbers"""
        scored = [self.scoreZmw(self.basH5[zmw]) for zmw in holeNumbers]
        return [self.makeLabeledZmw(scoreBunch) for scoreBunch in scored if scoreBunch.numAdapters]

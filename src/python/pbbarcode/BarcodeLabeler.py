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

def makeScoreAdaptersFunc(forwardScorer, reverseScorer, pairedScorer,
                          scoreMode, numSeqs, oldWorkflow):
    """Once the flanking regions around an adapter have been extracted
    they need to be scored against the appropriate set of barcode sequences,
    the specifics of which vary by scoreMode and workflow.
    """

    def scoreAdaptersOld(holeNum, adapters, scoredFirst):
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
                adapterScores[i] = np.maximum((fscores + rrscores)/1.0,
                                             (rscores + ffscores)/1.0)
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else np.zeros(numSeqs)

        return (holeNum, len(adapters), barcodeScores, adapterScores, scoredFirst)

    def scoreAdaptersPaired(holeNum, adapters, scoredFirst):
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
                adapterScores[i] = (fscores + rscores)/1.0
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else np.zeros(numSeqs)

        return (holeNum, len(adapters), barcodeScores, adapterScores, scoredFirst)

    def scoreAdaptersSymmetric(holeNum, adapters, scoredFirst):
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
                adapterScores[i] = (fscores + rscores)/1.0
            # Otherwise return the empty
            else:
                adapterScores[i] = barcodeScores

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else np.zeros(numSeqs)

        return (holeNum, len(adapters), barcodeScores, adapterScores, scoredFirst)

    # Return the selected scoreAdapters function
    if oldWorkflow:
        return scoreAdaptersOld
    elif scoreMode == 'paired' and not oldWorkflow:
        return scoreAdaptersPaired
    elif scoreMode == 'symmetric' and not oldWorkflow:
        return scoreAdaptersSymmetric

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
        self.barcodeLength   = np.unique(map(lambda x : len(x.sequence),
                                          self.barcodeFasta))
        self.useOldWorkflow  = useOldWorkflow
        self.adapterSidePad  = adapterSidePad
        self.insertSidePad   = insertSidePad
        self.maxHits         = maxHits
        self.scoreFirst      = scoreFirst
        self.startTimeCutoff = startTimeCutoff

        if scoreMode not in ['symmetric', 'paired']:
            raise Exception("scoreMode must either be symmetric or paired")
        self.scoreMode = scoreMode

        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

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

        # Make a "scoreAdapters" function for scoring flanking regions
        self.scoreAdapters = makeScoreAdaptersFunc(forwardScorer,
                                                   reverseScorer,
                                                   pairedScorer,
                                                   self.scoreMode,
                                                   self.numSeqs,
                                                   self.useOldWorkflow)

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
        adapters, scoredFirst = self._flankingSeqs(zmw)
        return self.scoreAdapters(zmw.holeNumber, adapters, scoredFirst)

    def labelZmws(self, holeNumbers):
        """Return a list of LabeledZmws for input holeNumbers"""

        def chooseSymmetric(o):
            """
            Convert a tuple "o" returned by scoreZmw into a LabeledZmw object
                for a symmetrically barcoded ZMW
            """
            p = np.argsort(-o[2])
            return LabeledZmw(o[0], o[1], p[0], o[2][p[0]], p[1], o[2][p[1]], o[3])

        def choosePaired(o):
            """
            Convert a tuple "o" returned by scoreZmw into a LabeledZmw object
                for an asymmetrically barcoded ZMW
            """
            if o[1] == 1:
                s = np.array([max(o[2][i], o[2][i + 1]) for i in \
                                 xrange(0, len(self.barcodeSeqs), 2)])
                p = np.argsort(-s)
                s = s[p]
            else:
                # score the pairs by scoring the two alternate
                # ways they could have been put on the molecule. A
                # missed adapter will confuse this computation.
                scores  = o[3]
                results = np.zeros(len(self.barcodeSeqs)/2)
                for i in xrange(0, len(self.barcodeSeqs), 2):
                    pths = [0,0]
                    for j in xrange(0, len(scores)):
                        pths[j % 2] += scores[j][i]
                        pths[1 - j % 2] += scores[j][i + 1]
                    results[i/2] = max(pths)

                p = np.argsort(-results)
                s = results[p]

            return LabeledZmw(o[0], o[1], p[0], s[0], p[1], s[1], o[3])

        # Select the "choose" method to be used formatting raw barcode scores
        #    into LabeledZmw Objects
        if self.scoreMode == 'symmetric':
            choose = chooseSymmetric
        elif self.scoreMode == 'paired':
            choose = choosePaired
        else:
            raise Exception("Unsupported scoring mode in BarcodeLabeler.py")

        scored = [self.scoreZmw(self.basH5[zmw]) for zmw in holeNumbers]

        return [choose(scoreTup) for scoreTup in scored if scoreTup[1]]

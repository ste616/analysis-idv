import json
import dominate
from dominate.tags import *

def startPage(title="Test page"):
    doc = dominate.document(title=title)

    with doc.head:
        link(rel='stylesheet', href='/people/ste616/ihv.css')

    with doc:
        div(title, _class='pagetitle')

    return doc

def endPage(doc, filename='test.html'):
    # Output the HTML file.
    with open(filename, 'w') as fp:
        print >>fp, doc

def outputIndex(epochData=None, outputDirectory=".", outputName=None,
                timescalePlots=None):
    if epochData is None:
        return
    if outputName is None:
        return

    # The page title will just be the source name.
    pageTitle = epochData[0]['timeTuple'][0]
    doc = startPage(pageTitle)

    # Build the page.
    epochsPerRow = 6
    with doc:
        with div(_class="epoch-container"):
            div("Available Epochs", _class="headingClass")
            with table():
                i = 0
                j = 0
                while i < len(epochData):
                    with tr():
                        j = 0
                        while j < epochsPerRow and i < len(epochData):
                            if 'timeTuple' in epochData[i]:
                                with td():
                                    a("%s (DOY %04d-%03d)" % (epochData[i]['epochName'],
                                                              epochData[i]['timeTuple'][2],
                                                              epochData[i]['timeTuple'][3]),
                                      href=epochData[i]['htmlTimescalesFile'])
                            i += 1
                            j += 1
            div("Timescale evolution plots", _class="headingClass")
            with div(_class="epoch-timescales"):
                with table():
                    with thead():
                        with tr():
                            td("Frequency")
                            td("Timescale variation")
                    with tbody():
                        for i in xrange(0, len(timescalePlots)):
                            with tr():
                                td("%d" % int(timescalePlots[i]['frequency']))
                                with td():
                                    with figure():
                                        img(src=timescalePlots[i]['plotFile'])

    # The filename to write to.
    pageFile = "%s/%s" % (outputDirectory, outputName)
    endPage(doc, pageFile)

def outputEpoch(epochData=None, outputDirectory=".", outputName=None, 
                linkNext=None, linkPrevious=None, linkUp=None):
    # Take some epoch of data, and make a web page that
    # displays some plots and some timescale quantities.
    if epochData is None or epochData['useable'] == False:
        return
    if outputName is None:
        return

    # The title will be the epoch name and source name.
    pageTitle = "%s (MJD %d DOY %04d-%03d)" % epochData['timeTuple']
    doc = startPage(pageTitle)

    # Build the page.
    with doc:
        with div(_class="epoch-container"):
            with div(_class="epoch-navigation"):
                if linkUp is not None:
                    a("Epoch Index", href=linkUp)
                if linkPrevious is not None:
                    a("Previous Epoch", href=linkPrevious)
                if linkNext is not None:
                    a("Next Epoch", href=linkNext)
            with div(_class="epoch-data"):
                if "dynamicImage" in epochData:
                    with figure():
                        img(src=epochData['dynamicImage'])
                        figcaption("Dynamic Spectrum")
                if "timeSeriesImage" in epochData:
                    with figure():
                        img(src=epochData['timeSeriesImage'])
                        figcaption("Time Series")
                if "timescaleFrequencyImage" in epochData:
                    with figure():
                        img(src=epochData['timescaleFrequencyImage'])
                        figcaption("Timescale vs Frequency")
            if 'acfPlotNames' in epochData:
                with div(_class="epoch-timescales"):
                    with table():
                        with thead():
                            with tr():
                                td("Frequency")
                                td("Auto-correlation function")
                                td("Timescale")
                                td("log10(Modulation Index)")
                        with tbody():
                            for i in xrange(0, len(epochData['modulationIndex']['frequencies'])):
                            #for i in xrange(0, len(epochData['acfPlotNames'])):
                                # Check if we have a plot to show.
                                showPlot = True
                                if i >= len(epochData['acfPlotNames']):
                                    showPlot = False
                                acfi = -1
                                if showPlot == True:
                                    for j in xrange(0, len(epochData['acfPlotNames'])):
                                        if epochData['acfPlotNames'][j]['frequency'] == int(epochData['modulationIndex']['frequencies'][i]):
                                            acfi = j
                                            break
                                with tr():
                                    if showPlot == True:
                                        #if epochData['acfPlotNames'][acfi]['frequency'] == 0:
                                        if acfi == -1:
                                            td("All")
                                        else:
                                            td("%d" % epochData['acfPlotNames'][acfi]['frequency'])
                                    else:
                                        td("%d" % int(epochData['modulationIndex']['frequencies'][i]))
                                    with td():
                                        if showPlot == True:
                                            with figure():
                                                img(src=epochData['acfPlotNames'][acfi]['fileName'])
                                    with td():
                                        if showPlot == True:
                                            #if epochData['acfPlotNames'][i]['frequency'] == 0:
                                            if acfi == -1:
                                                for j in xrange(0, len(epochData['acf']['timescale'])):
                                                    edat = epochData['acf']['timescale'][j]
                                                    if ((edat['value'] is not None) and 
                                                        (edat['valueError'] is not None) and
                                                        (edat['timeUnits'] is not None) and
                                                        (edat['timeUnits'] is not None)):
                                                        div("%.2f +/- %.2f %s at %d" % (edat['value'],
                                                                                        edat['valueError'],
                                                                                        edat['timeUnits'],
                                                                                        int(epochData['acf']['frequencies'][j])))
                                                    else:
                                                        div("NONE")
                                            else:
                                                #print epochData['acf']
                                                for j in xrange(0, len(epochData['acf']['frequencies'])):
                                                    edat = epochData['acf']['timescale'][j]
                                                    if (epochData['acf']['frequencies'][j] is not None 
                                                        and int(epochData['acf']['frequencies'][j]) == epochData['acfPlotNames'][i]['frequency']):
                                                        if (edat['value'] is None or
                                                            edat['valueError'] is None):
                                                            div("timescale undetermined at %d" % int(epochData['acf']['frequencies'][j]))
                                                        elif edat['lowerLimit'] == True:
                                                            div("> %.2f %s" % (edat['value'], edat['timeUnits']))
                                                        else:
                                                            div("%.2f +/- %.2f %s" % (edat['value'],
                                                                                    edat['valueError'],
                                                                                    edat['timeUnits']))
                                                            div("%d points, R=%.2f" % (edat['nFitPoints'], edat['fitQuality']))
                                                            if edat['degreesOfFreedom'] > 0:
                                                                div("fit chi-squared = %.2f (%.2f over %d)" % ((edat['chiSquared'] / float(edat['degreesOfFreedom'])), edat['chiSquared'], edat['degreesOfFreedom']))
                                                            else:
                                                                div("no points for chi-squared")
                                    with td():
                                        div("%.3f" % epochData['modulationIndex']['modulationIndex'][i])
                                                                               

    # The filename to write to.
    pageFile = "%s/%s" % (outputDirectory, outputName)

    endPage(doc, pageFile)


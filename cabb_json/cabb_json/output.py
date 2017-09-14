import json
import dominate
from dominate.tags import *

def startPage(title="Test page"):
    doc = dominate.document(title=title)

    with doc.head:
        link(rel='stylesheet', href='ihv.css')

    with doc:
        div(title, _class='pagetitle')

    return doc

def endPage(doc, filename='test.html'):
    # Output the HTML file.
    with open(filename, 'w') as fp:
        print >>fp, doc

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
                        with tbody():
                            for i in xrange(0, len(epochData['acfPlotNames'])):
                                with tr():
                                    if epochData['acfPlotNames'][i]['frequency'] == 0:
                                        td("All")
                                    else:
                                        td("%d" % epochData['acfPlotNames'][i]['frequency'])
                                    with td():
                                        with figure():
                                            img(src=epochData['acfPlotNames'][i]['fileName'])
                                    with td():
                                        if epochData['acfPlotNames'][i]['frequency'] == 0:
                                            for j in xrange(0, len(epochData['acf']['timescale'])):
                                                div("%.2f +/- %.2f %s at %d" % (epochData['acf']['timescale'][j]['value'],
                                                                                epochData['acf']['timescale'][j]['valueError'],
                                                                                epochData['acf']['timescale'][j]['timeUnits'],
                                                                                int(epochData['acf']['frequencies'][j])))
                                        else:
                                            #print epochData['acf']
                                            for j in xrange(0, len(epochData['acf']['frequencies'])):
                                                if (epochData['acf']['frequencies'][j] is not None 
                                                    and int(epochData['acf']['frequencies'][j]) == epochData['acfPlotNames'][i]['frequency']):
                                                    if (epochData['acf']['timescale'][j]['value'] is None or
                                                        epochData['acf']['timescale'][j]['valueError'] is None):
                                                        div("timescale undetermined at %d" % int(epochData['acf']['frequencies'][j]))
                                                    elif epochData['acf']['timescale'][j]['isRandom'] == True:
                                                        div("> %.2f %s" % (epochData['acf']['timescale'][j]['value'], epochData['acf']['timescale'][j]['timeUnits']))
                                                    else:
                                                        div("%.2f +/- %.2f %s" % (epochData['acf']['timescale'][j]['value'],
                                                                                  epochData['acf']['timescale'][j]['valueError'],
                                                                                  epochData['acf']['timescale'][j]['timeUnits']))
                                                        div("%d points, R=%.2f" % (epochData['acf']['timescale'][j]['nFitPoints'], epochData['acf']['timescale'][j]['fitQuality']))
                                                                               

    # The filename to write to.
    pageFile = "%s/%s" % (outputDirectory, outputName)

    endPage(doc, pageFile)


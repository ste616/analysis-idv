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
        fp.write(doc)

def outputEpoch(epochData=None, outputDirectory=".", linkNext=None, linkPrevious=None,
                linkUp=None):
    # Take some epoch of data, and make a web page that
    # displays some plots and some timescale quantities.
    if epochData is None or epochData['useable'] == False:
        return

    # The title will be the epoch name and source name.
    pageTitle = "%s (MJD %d DOY %04d-%03d)" % epochData['timeTuple']
    doc = startPage(pageTitle)

    # Build the page.
    with doc:
        with div(_class="epoch-container"):
            with div(_class="epoch-data"):
                if "dynamicImage" in epochData:
                    with figure():
                        img(src=epochData['dynamicImage'])
                        figcaption("Dynamic Spectrum")
                if "timeSeriesImage" in epochData:
                    with figure():
                        img(src=epochData['timeSeriesImage'])
                        figcaption("Time Series")


    # The filename to write to.
    pageFile = "%s/%s_doy%04d-%03d_epochinfo.html" % (outputDirectory,
                                                      epochData['timeTuple'][0],
                                                      epochData['timeTuple'][2],
                                                      epochData['timeTuple'][3])
    endPage(doc, pageFile)


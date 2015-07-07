"""
This module contains class that handles all data writing and saving to csv
file format for GenSIP.
"""
import csv
import GenSIP.functions as fun
import os
from socket import gethostname

###################################################################################

###################################################################################

class DataToCSV (object):
    
    def __init__(self,filepath,title,mode='w+b'):
        """ Creates a CSV file. 
                filepath - path of CSV file
                title - a string denoting the test being run. For cleantests, this
                    is the sample set string. 
            Kwargs:
                mode = 'w+b' - default set to allow writing."""
        assert type(title)==str, "Title must be a string"
        if filepath.endswith('.csv'):
            self.path = filepath
        elif os.path.isdir(filepath):
            self.path = os.path.join(filepath,title+'.csv')
        else:
            raise Exception("Path must be a directory or .csv file: {0}".format(filepath))
        
        self.CSVfile = open(self.path, mode)
        self.dataWriter = csv.writer(self.CSVfile)
        self.Title = title
        self.AllRows = []
        TitleRow = self.makeTitle(title)
        self.writeHeader(TitleRow)

    ###################################################################################
        
    def makeTitle(self, title):
        """
        Makes the title row (first row) of the CSV file based off of the title 
        vairable.
        Inputs: title string
        """
        TitleRow = ["GenSIP Data",'','','']
        # Special Case for sample set strings for cleaning tests
        if self.isCleaningTest(): 
                if title.endswith("NF"):
                    TitleRow[1] = "non-flight foils"
                    TitleRow[3] = title.strip("NF")+" Test"
                elif title.endswith("F"):
                    TitleRow[1] = "non-flight foils"
                    TitleRow[3] = title.strip("NF")+" Test"
        else:
            TitleRow[2] = title
            
        return TitleRow
        
    ###################################################################################

    def isCleaningTest(self):
        """
        Checks to see if the particular test is a cleaning test, by looking to see 
        if the title is all capitalized and ends with an F (as in BGF or BGNF
        """
        ret = self.Title.isupper() and self.Title.endswith("F")
        return ret
        
    ###################################################################################
    
    def writeHeader(self,TitleRow):
        """writes the header to the csv file given the first title row"""
        # Second row is the date and time of the run and the computer on which 
        # this function was run:
        datestring = fun.getDateString()
        host = os.path.splitext(gethostname())[0]
        Version = fun.getGenSIPVersion()
        InfoRow1 = ["Date:", datestring]
        InfoRow2 = ["Computer:", host]
        InfoRow3 = ["Last Modification:", Version]
        self.dataWriter.writerows([TitleRow,InfoRow1,InfoRow2,InfoRow3,['']])
        self.AllRows.extend([TitleRow,InfoRow1,InfoRow2,InfoRow3,['']])

    ###################################################################################

    def writeDataFromDict(self, dataDict, **kwargs):
        """
        Takes a dictionary argument with the data from an analysis fuction and 
        writes it to the csv file. Automatically creates column headers by 
        sorting alphabetically the individual pieces of data, or manually by 
        receiving a list describing the desired column headers. 
        This function should be easily generalized to receive any dictionary with the
        format:
            resultsDict = {
                        'item1':{'result1':<str or number>,
                                ...,
                                'resultN':<str or number>
                                },
                                
                                ...
                                
                        'itemN':{'result1':<str or number>,
                                ...,
                                'resultN':<str or number>
                                }  
                        'TOTALS':{'tot1':<entry1>,...,'totN':<entryN>}
                        'OTHER_FOOTER_INFO':{'other1:<entry1>,...,'otherN':<entryN>}
                        }   
                            
        In general 'item1' to 'itemN' will be the subimage names. 
        Inputs:
            - dataDict - results dictionary formatted as described above.
        Key-Word Arguments:
            - FirstColHead = 'Foil' - Header for the rows of items. 
            - footerItems = ["TOTALS","TOTAL"] - list of entries (keys in the dataDict)
                that are to only be included in the footer. Usually this entry is a 
                subDictionary of totals or averages to be placed at the end of the 
                csv file. 
            - colHeads = 0 - option to give the column headers of the data. If
                left at 0, writeDataFromDict will automatically create headers 
                from the result keys for each item. 
        """
        FirstColHead = kwargs.get('FirstColHead','Foil')
        notItems = kwargs.get('footerItems',["TOTALS","TOTAL"])
        colHeads = kwargs.get('colHeads',0)
        colHeadsProvided = type(colHeads) in [list, tuple, set]
        if not colHeadsProvided:
            colHeads = [FirstColHead]
        names = [k for k in dataDict.keys() if not(k in notItems)]
        names.sort()
        Rows = []
        # if the column headers are not given, generate them from the items in the
        # Data dictionary.
        if not colHeadsProvided:
            for sub in names:
                if type(dataDict[sub])==dict:
                    Keys = dataDict[sub].keys()
                    Keys.sort()
                    # Make sure all the keys in the subImg dictionary are 
                    # accounted for in the colHeads
                    colHeads.extend([k for k in Keys if (type(dataDict[sub][k])!=dict) and not(k in colHeads)])
        Rows.append(colHeads)
        """ Populate the rows for each item in the dictionary """
        for sub in names:
            if type(dataDict[sub])==dict:
                row = ["'"+sub]
                for col in colHeads[1:]:
                    entry = dataDict[sub].get(col,"'--")
                    row.append(entry)
            Rows.append(row)
            
        """Make the Footer"""
        # Separate the footer with a spacer row:
        spaceRow = ['']
        Rows.append(spaceRow)
        # Make footer by iterating through the footer items 
        footerItems = [k for k in notItems if k in dataDict.keys()]
        for cat in footerItems:
            row = [str(cat)+':']
            FooterRows = [row]
            if type(dataDict[cat])==dict:
                for k in dataDict[cat]:
                    footerEntry = [k,dataDict[cat][k]]
                    FooterRows.append(footerEntry)
            else:
                footerEntry = [dataDict[cat]]
                FooterRows.append(footerEntry)
            Rows.extend(FooterRows)
        self.dataWriter.writerows(Rows)
        self.AllRows.extend(Rows)
    
    ###################################################################################
    
    def closeCSVFile(self):
        """Closes the CSV file"""
        self.CSVfile.close()
        
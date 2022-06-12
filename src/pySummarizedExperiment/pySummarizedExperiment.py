import pandas as pd
import copy
import numpy as np

class pySummarizedExperiment:
    """
    A class for containing -omics style data without duplication
    
    Methods
    -------
    toLongDataFrame()
        Exports the pySummarizedExperiment to a single long-format DataFrame
    bind(exp, axis = 0)
        Bind two pySummarizedExperiments together. axis = 0 binds rowWise, axis = 1 binds columnWise 
    setIndex(index)
        Set a new rowIndex of the pySummarizedExperiment 
    index()
        Get the rowIndex of the pySummarizedExperiment
    setColnames(names)
        Set new column names of the pySummarizedExperiment 
    colnames()
        Get column names of the pySummarizedExperiment 
    assays()
        Get the assayData of the pySummarizedExperiment
    assay(i = 0)
        Get a specific assay, defaults to the first index of the keys
    colData(col = None)
        Get all the colData, or a specific column when given
    rowData(col = None)
        Get all the rowData, or a specific column when given
    metaData(key = None)
        Get all the metaData, or a specific value when a key is given 
    """
    def __init__(self, assays: dict = dict(), colData: pd.DataFrame = pd.DataFrame(), 
                 rowData: pd.DataFrame = pd.DataFrame(), metaData: dict = dict(),
                 longDf: pd.DataFrame = pd.DataFrame(), colIndex: str = None, rowIndex: str = None):
        """
        Initialize the object. This can be done in two ways:
        1) Supply assays, colData, rowData, and optional metaData.
        2) Supply a longDf, colIndex, and rowIndex
        
        Parameters
        ----------
        assays: dict
            A dictionary with Pandas DataFrames as values. Keys will be used as names for the assays
        colData: pd.DataFrame
            DataFrame typically used for samples. These samples should be the index of the DataFrame
        rowData: pd.DataFrame
            DataFrame typically used for features/genes. These features/genes should be the index of the DataFrame
        metaData: dict
            A dictionary with various values describing the experiment. Won't be altered by subsetting
        longDf: pd.DataFrame 
            DataFrame in the long-format with assays, samples & features. Can be used instead of assays, colData & rowData
        colIndex: str
            Which column is used as colIndex when using a long DataFrame
        rowIndex: str
            Which column is used as rowIndex when using a long DataFrame
        """

        if len(longDf.index) > 0:
            rowData, colData, assays, metaData = self.__longToExperiment(longDf, rowIndex, colIndex)
            
        self.__assayData = {str(key): assays[key] for key in assays.keys()}
        self.__colDat = colData
        self.__rowDat = rowData
        self.__metadat = metaData
        
    def toLongDataFrame(self):
        """
        Export the object to a Pandas DataFrame in long-format. Can be used to write to a file
        """
        longDf = pd.concat([x.stack().to_frame() for x in self.assays().values()], axis = 1)
        longDf.columns = self.assays().keys()
        longDf = pd.merge(longDf, self.rowData(), right_index=True, left_on = self.index().name)
        longDf = pd.merge(longDf, self.colData(), right_index=True, left_on = self.colData().index.name)
        return longDf.reset_index()
        
    def __longToExperiment(self, df, rowIndex, colIndex):
        metadata = np.unique(df[df.columns[np.where(df.nunique() == 1)]])
        df = df[df.columns[np.where(df.nunique() > 1)]]
        df[rowIndex] = df[rowIndex].astype(str)
        df[colIndex] = df[colIndex].astype(str)
        rowCols = [rowIndex]
        for col in df.columns:
            if col != rowIndex:
                comb = df[[rowIndex, col]].drop_duplicates()
                one = len(set(comb[rowIndex]))                
                if one == len(comb.index):
                    rowCols.append(col)
        rowData = df[rowCols].drop_duplicates().set_index(rowIndex)
        
                    
        colCols = [colIndex]
        for col in df.columns:
            if col != colIndex:
                comb = df[[colIndex, col]].drop_duplicates()
                one = len(set(comb[colIndex]))
                if one == len(comb.index):
                    colCols.append(col)
        colData = df[colCols].drop_duplicates().set_index(colIndex)

        assays = {name: df[[name, rowIndex, colIndex]].pivot(index = rowIndex, columns = colIndex, values = name) 
                  for name in df.columns if name not in rowCols + colCols}
        return rowData, colData, assays, metadata 
            
    def __validateWith(self, exp: "pySummarizedExperiment", axis: int):
        if axis == 0:
            cols = all([ind in exp.colnames() for ind in self.colnames()])
            cols2 = all([ind in exp.colData().columns for ind in self.colData().columns])

            rows = all([ind not in exp.index() for ind in self.index()])
            return cols and cols2 and rows
        
        cols = all([ind not in exp.colnames() for ind in self.colnames()])
        rows = all([ind in exp.index() for ind in self.index()])
        rows2 = all([ind in exp.rowData().columns for ind in self.rowData().columns])

        return cols and rows and rows2
    
    def __bindRowData(self, exp: "pySummarizedExperiment"):
        self.__rowDat = pd.concat([self.rowData(), exp.rowData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 0)        

    def __bindColData(self, exp: "pySummarizedExperiment"):
        self.__colDat =  pd.concat([self.colData(), exp.colData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 1)   
            
    def __inPlace(self, inPlace):
        if not inPlace:
            x = copy.deepcopy(self)
        else: 
            x = self
        return x    
         
        
    def bind(self, exp: "pySummarizedExperiment" , axis: int = 0, inPlace: bool = False):
        """
        Bind another pySummarizedExperiment
        
        Parameters
        ----------
        exp: pySummarizedExperiment
            Initialized pySummarizedExperiment with either the same rows or same columns depending on the axis. 
        axis: int
            Which axis should be used for binding? Either 0 (rowWise), or 1 (columnWise)
        inPlace: bool
            Should the binding be in-place (True) or returning a copy (False, default) 
        """
        if axis not in [0, 1]:
            raise Exception("axis must be in range [0, 1]")
        
        x = self.__inPlace(inPlace) 
            
        if x.__validateWith(exp, axis):
            if axis == 0:
                x.__bindRowData(exp)
            else: 
                x.__bindColData(exp)
        return x

        
    def __rownames(self, i: str):
        self.__assayData = {key: x.loc[i] for key, x in self.__assayData.items()}
        self.__rowDat = self.rowData().loc[i]
        return self

    def __columns(self, i: str):      
        self.__assayData = {key: x[i] for key, x in self.__assayData.items()}
        self.__colDat = self.colData().loc[i]
        return self

    def setIndex(self, index: str, inPlace: bool = False):
        """
        Set the row index of the pySummarizedExperiment
        
        Parameters
        ----------
        
        index: str
            Row index to be set 
        inPlace: bool
            Should the replacement be in-place (True) or returning a copy (False, default)         
        """
        x = self.__inPlace(inPlace) 
        
        if len(index) == len(x.index()):
            x.__rowDat.index = index
            for key in x.__assayData.keys():
                x.__assayData[key].index = index
        return x
        

    def index(self):
        """
        Get the row index of the pySummarizedExperiment
        """
        return self.rowData().index
    
    def setColnames(self, names: str, inPlace: bool = False):
        """
        Set the column index of the pySummarizedExperiment
        
        Parameters
        ----------
        names: str
            Names of the new columns to use.
        inPlace: bool
            Should the replacement be in-place (True) or returning a copy (False, default)         

        """
        x = self.__inPlace(inPlace) 
        
        if len(names) == len(x.colnames()):
            x.__colDat.index = names
            for key in x.__assayData.keys():
                x.__assayData[key].columns = names
        return x

    
    def colnames(self):
        """
        Get the column index (names) of the pySummarizedExperiment
        """
        return self.colData().index

    def __setitem__(self, i, assay: pd.DataFrame):
        # check if assay dimnames == assayData dimnames
        self.assayData[i] = assay
    
    def __parseArgs(self, arg, lookup):
        if type(arg) == pd.core.series.Series:
            arg = np.where(arg)[0].tolist()
        if type(arg) == slice:
            arg = lookup[arg]
        elif type(arg) != list:
            if type(arg) == int:
                arg = lookup[arg]
            arg = [arg]
        elif type(arg) == list:
            for j in range(len(arg)):
                if type(arg[j]) == int:
                    arg[j] = lookup[arg[j]]
        return arg

    def __getitem__(self, i):
        try:
            rows, cols = i        
        except ValueError:
            rows = slice(None)
            cols = i
        
        rows = self.__parseArgs(rows, self.index())
        cols = self.__parseArgs(cols, self.colnames())
                
        
        if type(rows) != list and type(rows) != pd.core.indexes.base.Index:
            rows = [rows]
        if type(cols) != list and type(cols) != pd.core.indexes.base.Index:
            cols = [cols]    
        
        return copy.copy(self).__rownames(rows).__columns(cols)

    def assays(self):
        """
        Get the dictionary with stored assays 
        """
        return self.__assayData

    def assay(self, i = 0):
        """
        Get a specific assay by key or index
        
        Parameters
        ----------
        i: int, str
            Key (str) or index (int) of the assay to retrieve
        """
        if type(i) == int:
            i = list(self.__assayData.keys())[i]
        return self.__assayData[i]
   
    def colData(self, col: str = None):
        """
        Retrieve all colData or when given, a specific column
        
        Parameters
        ----------
        col: str
            Optional column to retrieve from the rowData
        """
        if col == None:
            col = self.__colDat.columns
        return self.__colDat[col]
    
    def rowData(self, col: str = None):
        """
        Retrieve all rowData or when given, a specific column
        
        Parameters
        ----------
        col: str
            Optional column to retrieve from the rowData
        """
        if col == None:
            col = self.__rowDat.columns
        return self.__rowDat[col]

    def metaData(self, key = None):
        """
        Retrieve all metaData or when given, a specific value
        
        Parameters
        ----------
        key: str
            Optional value to retrieve from the metaData
        """            
        if key == None:
            key = self.__metadat.keys()
        return self.__metadat[key]
    
    def __repr__(self) -> str:
        s = ["class: pySummarizedExperiment",
        "dim: %s %s" % (len(self.index()), len(self.colnames())),
        "metadata(%d): %s" % (len(self.metaData()), self.metaData().keys()), 
        "rownames(%d): %s" % (len(self.index()), ", ".join(self.index()[0:5].astype(str))),
        "rowData names(%d): %s" % (len(self.rowData().columns),  ", ".join(self.rowData()[0:5].astype(str))),
        "colnames(%d): %s" % (len(self.colnames()), ", ".join(self.colnames()[0:5].astype(str))),
        "colData names(%d): %s" % (len(self.colData().columns), ", ".join(self.colData()[0:5].astype(str))),
        "assays(%d): %s" % (len(self.assays()), ", ".join(list(self.assays().keys())[0:5])) 
        ]
        return "\n".join(s)
    



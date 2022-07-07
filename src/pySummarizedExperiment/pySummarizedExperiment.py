import sys
from matplotlib.pyplot import autoscale
import pandas as pd
import copy
import numpy as np
from datetime import datetime
from sklearn.decomposition import PCA

class pySummarizedExperiment:
    """
    A class for containing -omics style data without duplication. Based on `pandas` syntax.
    
    ## Parameters
    Initializing the object can be done in two ways:
    
    1. Supply assays, colData, rowData, and optional metaData.
    2. Supply a longDf, colIndex, rowIndex and optional metaData
        
    ### First option
    assays (`dict`): A dictionary with one ore more DataFrame objects as values. Keys will be used as names for the assays
    colData (`pd.DataFrame`): DataFrame typically used for samples / cells / patients. It is recommend to use identifiers of these as indexes of the DataFrame
    rowData (`pd.DataFrame`): DataFrame typically used for compound / genes. It is recommend to use identifiers of these as indexes of the DataFrame 
    metaData (`dict`): (optional) Dictionary with various values describing the experiment. Won't be altered by subsetting
    
    ### Second option
    longDf (`pd.DataFrame`): DataFrame in the long-format with assays, samples & features. Can be used instead of assays, colData & rowData, but should be supplied with `colIndex` and `rowIndex`.
    colIndex (`str`): Name of the column to be used from `longDf`. These values will be used as columns in the assay and as index in the colData.
    rowIndex (`str`): Name of the column to be used from `longDf`. These values will be used as columns in the assay and as index in the rowData.
    
    """
    def __init__(self, assays: dict = dict(), colData: pd.DataFrame = pd.DataFrame(), 
                 rowData: pd.DataFrame = pd.DataFrame(), metaData: dict = dict(),
                 longDf: pd.DataFrame = pd.DataFrame(), colIndex: str = None, rowIndex: str = None):

        if len(longDf.index) > 0:
            longDf = longDf[longDf.columns[np.where(longDf.nunique() > 1)]]
            rowData, colData, assays = self.__longToExperiment(longDf, rowIndex, colIndex)
            
            
        self.__assayData = {str(key): assays[key] for key in assays.keys()}
        self.__colDat = colData
        self.__rowDat = rowData
        self.__metadat = metaData
        self.__metadat["Datetime"] = datetime.now()
        
    def toLongDataFrame(self):
        """
        Export the object to a Pandas DataFrame in long-format. Can be used to write to a file
        """
        longDf = pd.concat([x.stack().to_frame() for x in self.assays().values()], axis = 1)
        longDf.columns = self.assays().keys()
        longDf = pd.merge(longDf, self.rowData(), right_index=True, left_on = self.index().name)
        longDf = pd.merge(longDf, self.colData(), right_index=True, left_on = self.colData().index.name)
        return longDf.reset_index()
    
    def __extractFromDf(self, df, index):
        rowCols = [index]
        for col in df.columns:
            if col != index:
                comb = df[[index, col]].drop_duplicates()          
                if len(set(comb[index])) == len(comb.index):
                    rowCols.append(col)
        return df.loc[:, rowCols].drop_duplicates().set_index(index)
        
    def __longToExperiment(self, df, rowIndex, colIndex):
        df[rowIndex] = df.loc[:, rowIndex].astype(str)
        df[colIndex] = df.loc[:, colIndex].astype(str)
        
        rowData = self.__extractFromDf(df, rowIndex)
        colData = self.__extractFromDf(df, colIndex)

        used = list(rowData.columns) + list(colData.columns) + [rowIndex, colIndex]
        assays = {name: df[[name, rowIndex, colIndex]].pivot(index = rowIndex, columns = colIndex, values = name) 
                  for name in df.columns if name not in used}
        return rowData, colData, assays 
            
    def __validateWith(self, exp: "pySummarizedExperiment", axis: int) -> bool:
        if axis == 0:
            cols = all([ind in exp.colnames() for ind in self.columns()])
            cols2 = all([ind in exp.colData().columns for ind in self.colData().columns])

            rows = all([ind not in exp.index() for ind in self.index()])
            return cols and cols2 and rows
        
        cols = all([ind not in exp.colnames() for ind in self.columns()])
        rows = all([ind in exp.index() for ind in self.index()])
        rows2 = all([ind in exp.rowData().columns for ind in self.rowData().columns])

        return cols and rows and rows2
    
    def __bindRowData(self, exp: "pySummarizedExperiment") -> None:
        self.__rowDat = pd.concat([self.rowData(), exp.rowData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 0)        

    def __bindColData(self, exp: "pySummarizedExperiment") -> None:
        self.__colDat =  pd.concat([self.colData(), exp.colData()])
        for key in self.__assayData.keys():
            self.__assayData[key] = pd.concat([self.assay(key), exp.assay(key)], axis = 1)   
         
        
    def bind(self, exp: "pySummarizedExperiment" , axis: int = 0) -> "pySummarizedExperiment":
        """
        Bind another pySummarizedExperiment
        
        Parameters:
            exp (pySummarizedExperiment): Initialized pySummarizedExperiment with either the same rows or same columns depending on the axis. 
            axis (int): Which axis should be used for binding? Either 0 (rowWise), or 1 (columnWise)
        """
        if axis not in [0, 1]:
            raise Exception("axis must be in range [0, 1]")
        
            
        if self.__validateWith(exp, axis):
            if axis == 0:
                self.__bindRowData(exp)
            else: 
                self.__bindColData(exp)
        return self

        
    def __rownames(self, i: str) -> "pySummarizedExperiment":
        self.__assayData = {key: x.loc[i] for key, x in self.__assayData.items()}
        self.__rowDat = self.rowData().loc[i]
        return self

    def __colnames(self, i: str) -> "pySummarizedExperiment":      
        self.__assayData = {key: x[i] for key, x in self.__assayData.items()}
        self.__colDat = self.colData().loc[i]
        return self

    def setIndex(self, index: list) -> "pySummarizedExperiment":
        """
        Set the row index of the pySummarizedExperiment
        
        Parameters:        
            index (list(str)): Row index to be set         
        """
        
        if len(index) == len(self.index()) and len(index) == len(set(index)):
            self.__rowDat.index = index
            for key in self.__assayData.keys():
                self.__assayData[key].index = index
        return self
    
    def setColumns(self, columns: list) -> "pySummarizedExperiment":
        """
        Set the columns and colData index of the pySummarizedExperiment
        
        Parameters:
            columns(list(str)): columns to be set      
        """
        if len(columns) == len(self.columns()) and len(columns) == len(set(columns)):
            self.__colDat.index = columns
            for key in self.__assayData.keys():
                self.__assayData[key].columns = columns
        return self
    
    def columns(self) -> pd.Index:
        """
        Get the column index (names) of the pySummarizedExperiment
        """
        return self.colData().index
        

    def index(self) -> pd.Index:
        """
        Get the row index of the pySummarizedExperiment
        """
        return self.rowData().index
    
    def log2(self, assay = 0, df: pd.DataFrame = None) -> pd.DataFrame:
        if df is None:
            df = self.assay(assay)
        return pd.DataFrame(np.log2(df))
    
    def autoScale(self, assay = 0, df: pd.DataFrame = None) -> pd.DataFrame:
        """
        Auto-Scale an assay
        
        Parameters:
            assay (int / str): Index or Name of assay to be auto-scaled
        """
        if df is None:
            df = self.assay(assay)
        return pd.DataFrame((df - df.mean()) / df.std())
    
    def __preprocessing(self, assay, log2, autoscale, transpose):
        if log2:
            df = self.log2(assay)
        else:
            df = self.assay(assay)

        if autoscale:
            df = self.autoScale(df = df)
        
        if transpose:
            df = df.transpose()
            
        return df
    
    def boxplot(self, assay = 0, log2 = False, autoscale = False, transpose = False):
        return self.__preprocessing(assay, log2, autoscale, transpose).boxplot()
        
    
    def pca(self, assay = 0, components = 5, log2 = False, autoscale = False, transpose = True) -> pd.DataFrame:
        df = self.__preprocessing(assay, log2, autoscale, transpose)
        pca = PCA(n_components=components)
        df = pca.fit_transform(df)
        return pd.DataFrame(df)

    def __validate_assay(self, assay: pd.DataFrame) -> bool:
        return all([x in assay.index for x in self.index()]) and \
            all([x in assay.columns for x in self.columns()])

    def __setitem__(self, i, assay: pd.DataFrame) -> None:
        if self.__validate_assay(assay):
            self.__assayData[i] = assay
    
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

    def __getitem__(self, i) -> "pySummarizedExperiment":
        """
        Subset the pySummarizedExperiment object by slices, masks or names
        
        Parameters:
            i (list / slice / mask): Indexes to subset
        """
        try:
            rows, cols = i        
        except ValueError:
            rows = slice(None)
            cols = i
        
        
        rows = self.__parseArgs(rows, self.index())
        cols = self.__parseArgs(cols, self.columns())
                
        
        if type(rows) != list and type(rows) != pd.core.indexes.base.Index:
            rows = [rows]
        if type(cols) != list and type(cols) != pd.core.indexes.base.Index:
            cols = [cols]    
        
        return copy.copy(self).__rownames(rows).__colnames(cols)

    def assays(self) -> dict:
        """
        Get the dictionary with stored assays 
        """
        return self.__assayData

    def assay(self, i = 0) -> pd.DataFrame:
        """
        Get a specific assay by key or index
        
        Parameters:
            i (int / str): Key (str) or index (int) of the assay to retrieve
        """
        if type(i) == int:
            i = list(self.__assayData.keys())[i]
        return self.__assayData[i]
   
    def colData(self, column: str = None) -> pd.DataFrame:
        """
        Retrieve all colData or when given, a specific column
        
        Parameters:
            column (str): Optional column to retrieve from the rowData
        """
        if column == None:
            column = self.__colDat.columns
        return self.__colDat.loc[:, column]
    
    def rowData(self, column: str = None) -> pd.DataFrame:
        """
        Retrieve all rowData or when given, a specific column
        
        Parameters:
            column (str): Optional column to retrieve from the rowData
        """
        if column == None:
            column = self.__rowDat.columns
        return self.__rowDat.loc[:, column]

    def metaData(self, key = None):
        """
        Retrieve all metaData or when given, a specific value
        
        Parameters:
            key (str): Optional value to retrieve from the metaData
        """            
        if key is None:
            return self.__metadat
        return self.__metadat[key]
    
    
        
    
    def __sizeof__(self) -> int:
        assays = sum([sys.getsizeof(value) for key, value in self.assays().items()])
        return sys.getsizeof(self.colData()) + sys.getsizeof(self.rowData()) + assays
    
    def __repr__(self) -> str:
        s = ["class: pySummarizedExperiment",
        "dim: %s %s" % (len(self.index()), len(self.columns())),
        "metadata(%d): %s" % (len(self.metaData().keys()), ", ".join(list(self.metaData().keys())[0:5])), 
        "rownames(%d): %s" % (len(self.index()), ", ".join(self.index()[0:5].astype(str))),
        "rowData names(%d): %s" % (len(self.rowData().columns),  ", ".join(self.rowData()[0:5].astype(str))),
        "colnames(%d): %s" % (len(self.columns()), ", ".join(self.columns()[0:5].astype(str))),
        "colData names(%d): %s" % (len(self.colData().columns), ", ".join(self.colData()[0:5].astype(str))),
        "assays(%d): %s" % (len(self.assays()), ", ".join(list(self.assays().keys())[0:5])) 
        ]
        return "\n".join(s)
    



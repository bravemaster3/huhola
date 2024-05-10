#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script defines a HuHoLa class for performing a series of operations on digital elevation models (DEMs)
to generate microtopography classifications.

Usage:
    - Initialize the class with necessary parameters.
    - Use the `generate_microtopography` method to execute the entire workflow.

Example:
    from microtopography import HuHoLa #This might need appending the huhola package to path using sys
    m = HuHoLa(path_wbt = "path/to/whiteboxtools/binaries",
           wd = "path/to/working_dir",
           path_dem = "path/to/the_DEM.tif",
           threshold_fill = 0.0)
    m.generate_microtopography()
    #Then, you can change working directory if needed, and classify with different thresholds using the classify_microtopography method
    m.wd = "path/to/new_working_dir"
    m.threshold_fill = 0.04
    m.classify_microtopography()

Outputs:
    - DEM_filled.tif: filled DEM (aka hollow layer)
    - DEM_inverted_filled.tif: filled inverted DEM (aka hummock layer)
    - hol_hum_depth_height.tif: filled height (negative values) and depths (positive values) of hummocks and hollows
    - FINAL_MICROTOPO3.tif: 3 classes microtopography (lawn: 0, Hollow: 1, Hummock: 2)
    - FINAL_MICROTOPO5.tif: 5 classes microtopography (lawn: 0, Hollow: 1, Hummock: 2, lower-level lawns: 3, upper-level lawns: 4)
    
Created on Tue Nov 21 08:24:28 2023

@authors: Koffi Dodji Noumonvi <koffi.noumonvi@slu.se> & Nils Helge Havertz <nils.havertz@stud.uni-greifswald.de>
"""

import os
import whitebox_tools as wb
import rasterio
import numpy as np
from rasterio.features import sieve
import geopandas as gpd

import csv
import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import os


import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.metrics import cohen_kappa_score, confusion_matrix, accuracy_score, mean_squared_error, classification_report




class InvalidPath(Exception):
    """InvalidPath Exception class"""""
    pass

class HuHoLa():
    """
    Main class for the HuHoLa analysis. It performs operations such as DEM filling, microtopography classification, and sieving.
    
    Parameters:
        path_wbt (str): Path to the Whitebox Tools directory.
        wd (str): Working directory.
        path_dem (str): Path to the digital elevation model (DEM) file.
        fill_method (str, optional): Method for filling depressions in the DEM. Default is "wang_liu".
        fix_flats (bool, optional): Whether to fix flats during DEM filling. Default is False.
        fix_flats_inv (bool, optional): Whether to fix flats during the filling of the inverted DEM. Default is False.
        flat_increment (float, optional): Increment value for flat areas during DEM filling. Default is 0.0001.
        threshold_fill (float, optional): Threshold value for classifying microtopography features. Default is 0.0.
        number_classes (int, optional): Number of classes for microtopography classification. Default is None.
        sieve (bool, optional): Whether to perform sieving on the final classification results. Default is False.
        min_sieve_size (int, optional): Minimum number of pixels to retain during sieving. Any pixel group below the threshold is merged with its surrounding  Default is 2.
    """
    def __init__(self,
                 path_wbt: str,
                 wd: str,
                 path_dem: str,
                 fill_method: str = "wang_liu",
                 fix_flats: bool = False,
                 fix_flats_inv: bool = False,
                 flat_increment: float = 0.001,
                 threshold_fill: float = 0.0,
                 number_classes: int = None,
                 sieve: bool = False,
                 min_sieve_size: int = 2):
        """Instanciation of the class"""
        
        # checking that pp_wbt is a valid dir and 
        # create a white boxtools instance
        wbt = wb.WhiteboxTools()
        
        # Setting path_wbt using its setter method
        self.path_wbt = path_wbt
        wbt.set_whitebox_dir(self.path_wbt)
        
        # setting wd using its setter method
        self.wd = wd
        
        # setting the path_dem using its setter method
        self.path_dem = path_dem
        
        # Creating a dict to map the fill methods to the right methods
        self.__fill_method_dict = {'wang_liu': wbt.fill_depressions_wang_and_liu,
                            'simple': wbt.fill_depressions, 
                            'planchon_and_darboux': wbt.fill_depressions_planchon_and_darboux}
        
        # Setting the fill method protected attribute
        self.fill_method = fill_method
        
        # Setting threshold_fill using the setter method
        self.threshold_fill = threshold_fill
        
        
        # Setting the number of classes protected attribute. requires self._threshold_fill to be set  
        self.number_classes = number_classes
        
        # Setting the fix flats and flat increment using the setter methods
        self.fix_flats = fix_flats
        self.fix_flats_inv = fix_flats_inv
        self.flat_increment = flat_increment
        
        # Setting the sieve and min_sieve_size using their setter methods
        self.sieve = sieve
        self.min_sieve_size = min_sieve_size
        
    
    @property
    def path_wbt(self):
        """Getter for path to whitebox tools protected attribute _path_wbt"""
        return self._path_wbt
    
    @path_wbt.setter
    def path_wbt(self, value):
        """Setter for path to whitebox tools protected attribute _path_wbt"""
        if type(value) is str and os.path.isdir(value):
            try:
                self._path_wbt = value
            except Exception as err:
                raise InvalidPath(err)
        else:
            raise InvalidPath(f"path to Whitebox tools binaries '{value}' must be a valid directory")
    
    @property
    def wd(self):
        """Getter for working directory protected attribute _wd"""
        return self._wd
    
    @wd.setter
    def wd(self, value):
        """Setter for working directory protected attribute _wd"""
        if type(value) is str and os.path.isdir(value):
            self._wd = value
            #Creating output file paths
            self._p_dem_prep = os.path.join(self.wd, "DEM_prep.tif")
            self._p_filled_dem = os.path.join(self.wd, "DEM_filled.tif")
            self._p_hollow_lyr = os.path.join(self.wd, "DEM_subtr.tif")
            self._p_inv_dem = os.path.join(self.wd, "DEM_inverted.tif")
            self._p_filled_inv_dem = os.path.join(self.wd, "DEM_inverted_filled.tif")
            self._p_hummock_lyr = os.path.join(self.wd, "DEM_inverted_filled_subtr.tif")
            self._p_hol_hum_depth_height = os.path.join(self.wd, "hol_hum_depth_height.tif")
            
            # 3 class microtopography
            self._p_final_microtopo3 = os.path.join(self.wd, "FINAL_MICROTOPO3.tif")
            # 5 class microtopography
            self._p_final_microtopo5 = os.path.join(self.wd, "FINAL_MICROTOPO5.tif")
            
            # 3 class microtopography
            self._p_final_microtopo3_sieved = os.path.join(self.wd, "FINAL_MICROTOPO3_SIEVED.tif")
            # 5 class microtopography
            self._p_final_microtopo5_sieved = os.path.join(self.wd, "FINAL_MICROTOP5_SIEVED.tif")
        else:
            raise InvalidPath(f"working dir '{value}' must be a valid directory")
    
    @property
    def path_dem(self):
        """Getter for the digital elevation model path attribue"""
        return self._path_dem
    
    @path_dem.setter
    def path_dem(self, value):
        """Setter for the digital elevation model path attribute"""
        if type(value) is str and os.path.isfile(value):
            self._path_dem = value
        else:
            raise InvalidPath(f"DEM '{value}' must be a valid file")
    
    @property
    def fill_method(self):
        """Getter for the fill method attribute"""
        return self._fill_method
    
    @fill_method.setter
    def fill_method(self, value):
        if value in list(self.__fill_method_dict.keys()):
            self._fill_method = value
        else:
            print("Only three fill methods are supported: "
                  "Please choose 'wang_liu', 'simple', or 'planchon_and_darboux' "
                  "defaulting to the default method: Wang Liu")
            self._fill_method = "Wang_Liu"
    
    @property
    def threshold_fill(self):
        """Getter for threshold_fill"""
        return self._threshold_fill
    
    @threshold_fill.setter
    def threshold_fill(self, value):
        """Setter for threshold_fill"""
        if type(value) is not float:
            raise TypeError("threshold_fill must be a float [recommended values between 0.00 and 0.07, higher with better DEM resolutions]")
        elif value < 0.0:
            raise ValueError("threshold_fill  must be >= 0.0")
        elif value > 0.07:
            print("!!Warning!!: High threshold fill. Double check that this is desired." 
                  "Ignore this warning. The value will be set")
            self._threshold_fill = value
        else:
            self._threshold_fill = value
    
    @property
    def number_classes(self):
        """Getter method for protected attribute _number_classes"""
        return self._number_classes
    
    @number_classes.setter
    def number_classes(self, value):
        """Setter method for protected attribute _number_classes"""
        if value is None and self.threshold_fill == 0.0:
            print("changing the number_classes to 3 because threshold_fill is 0.0 and 3 classes and 5 classes are the same")
            self._number_classes = 3
        elif value in [3, 5, None]:
            self._number_classes = value
        else:
            raise ValueError("Only None => defaults to 3 & 5, 3 or 5 classes are supported")
            
    @property
    def fix_flats(self):
        """Getter for fix_flats"""
        return self._fix_flats
    
    @fix_flats.setter
    def fix_flats(self, value):
        """Setter for fix_flats"""
        if type(value) is not bool:
            raise TypeError("fix_flats must be True or False")
        else:
            self._fix_flats = value

    @property
    def fix_flats_inv(self):
        """Getter for fix_flats_inv"""
        return self._fix_flats_inv
    
    @fix_flats_inv.setter
    def fix_flats_inv(self, value):
        """Setter for fix_flats_inv"""
        if type(value) is not bool:
            raise TypeError("fix_flats_inv must be True or False")
        else:
            self._fix_flats_inv = value
    
    
    @property
    def flat_increment(self):
        """Getter for flat_increment"""
        return self._flat_increment
    
    @flat_increment.setter
    def flat_increment(self, value):
        """Setter for flat_increment"""
        if type(value) is not float:
            raise TypeError("flat_increment must be a float [recommended values between 0.0001 and 0.01, higher with better DEM resolutions]")
        elif value < 0.0:
            raise ValueError("flat_increment  must be >= 0.00")
        elif value > 0.01:
            print("!!Warning!!: High flat_increment. Double check that this is desired." 
                  "Ignore this warning. The value will be set")
            self._flat_increment = value
        else:
            self._flat_increment = value

    @property
    def sieve(self):
        """Getter for sieve"""
        return self._sieve
    
    @sieve.setter
    def sieve(self, value):
        """Setter for sieve"""
        if type(value) is not bool:
            raise TypeError("sieve must be True or False")
        else:
            self._sieve = value
    
    @property
    def min_sieve_size(self):
        """Getter method for protected attribute _min_sieve_size"""
        return self._min_sieve_size
    
    @min_sieve_size.setter
    def min_sieve_size(self, value):
        """Setter method for protected attribute _min_sieve_size"""
        if type(value) is not int:
            raise TypeError("min_sieve_size must be an integer")
        elif value < 0:
            raise ValueError("min_sieve_size must be >= 0")
        else:
            self._min_sieve_size = value
    
        
                
    def set_nodata_to_zero(self):
        """Sets the nodata to 0 in the DEM before proceeding with the rest of the analysis"""
        with rasterio.open(self.path_dem) as src:
            data = src.read(1, masked=True)
            nodata_value = src.nodata  # Get the nodata value from the original raster
            if nodata_value is not None:
                # Check if nodata_value is present in the data
                if np.ma.is_masked(data) and nodata_value in data.data:
                    print(f"Nodata value found in the original data: {nodata_value}")
                    # Set nodata to 0
                    data = np.ma.filled(data, fill_value=0)
                else:
                    print("No nodata values found in the original data.")
    
            meta = src.meta
            meta.update(nodata=0)
            with rasterio.open(self._p_dem_prep, 'w', **meta) as dst:
                dst.write(data, 1)

    
    def fill_dem(self):
        """Fills and saves the filled dem and filled inverted dem"""
        
        #Filling of the DEM
        print("Fill 1: ...DEM...")
        self.__fill_method_dict[self.fill_method](
            dem = self._p_dem_prep, 
            output = self._p_filled_dem, 
            fix_flats=self.fix_flats, 
            flat_increment=self.flat_increment)
        print("End of Fill 1: ...DEM...")
        
        #Let's read the DEM_prep again, saving its meta in self for future uses
        with rasterio.open(self._p_dem_prep, "r") as src:
            self.meta = src.meta
            dem = src.read()
            dem[dem == 0] = np.nan
        
        #Updating the meta
        self.meta.update({'nodata' : np.nan})
            
        with rasterio.open(self._p_filled_dem, "r") as src:
            dem_filled = src.read()
            dem_filled[dem_filled < -500] = np.nan #The non filled pixels were nodata, but we need them to be 0 instead
            dem_filled[dem_filled == 0] = np.nan
        
        src.close()
        
        #replacement of the dem_filled in the working directory as well
        #os.unlink(self._p_filled_dem)
        with rasterio.open(self._p_filled_dem,"w",**self.meta) as dest:
            dest.write(dem_filled)
        
        # Saving the filled depth of the hollows in the instance and writing also to disk
        self.hollows_lyr = dem_filled - dem
        with rasterio.open(self._p_hollow_lyr, "w", **self.meta) as dest:
            dest.write(self.hollows_lyr)

        #Inverting DEM by subtracting the dem from its maxvalue 
        dem_inv = np.nanmax(dem) - dem

        #Exporting the inverted DEM
        with rasterio.open(self._p_inv_dem, "w", **self.meta) as dest:
            dest.write(dem_inv)
            
        #*********************************************************#
        #Filling the inverted DEM
        print("Fill 2: ...Inverted DEM...")
        self.__fill_method_dict[self.fill_method](
            dem = self._p_inv_dem, 
            output = self._p_filled_inv_dem, 
            fix_flats = self.fix_flats_inv, 
            flat_increment = self.flat_increment)
        print("End of Fill 2: ...Inverted DEM...")
        
        with rasterio.open(self._p_filled_inv_dem, "r") as src:
            dem_inv_filled = src.read()
            dem_inv_filled[dem_inv_filled < -500] = np.nan #The non filled pixels were nodata, but we need them to be 0 instead
            dem_inv_filled[dem_inv_filled == 0] = np.nan
                
        src.close()
        
        #os.unlink(self._p_filled_inv_dem)
        with rasterio.open(self._p_filled_inv_dem,"w",**self.meta) as dest:
            dest.write(dem_inv_filled)

        self.hummock_lyr = dem_inv_filled - dem_inv

        
        with rasterio.open(self._p_hummock_lyr, "w", **self.meta) as dest:
            dest.write(self.hummock_lyr)
            
        #Direct difference =>  hol_hum_depth_height: positive = hollow, negative = hummocks 
        self.hol_hum_depth_height = self.hollows_lyr - self.hummock_lyr
        with rasterio.open(self._p_hol_hum_depth_height, "w", **self.meta) as dest:
            dest.write(self.hol_hum_depth_height)
    
    def classify_microtopography(self):
        """Uses the difference between hummock and hollow layer, and classification into hummock, hollow and lawn"""
        
        # #Before anything, let's refresh the init call because this method can be run separately and needs some variables to be updated
        # all_vars = self.__dict__.copy()
        # all_vars = {key: value for key, value in all_vars.items() if not key.startswith('p_') \
        #             and not key.endswith('_lyr') and key not in ['meta', 'fill_method_dict']}
        
        # self.__init__(**all_vars)
        # self.number_classes = None
        
        # #Direct difference approach
        # hol_hum_depth_height = self.hollows_lyr - self.hummock_lyr
        # with rasterio.open(self._p_hol_hum_depth_height, "w", **self.meta) as dest:
        #     dest.write(hol_hum_depth_height)
        
        hol_hum = self.hol_hum_depth_height.copy()
        hol_hum_depth_height = self.hol_hum_depth_height.copy()
        #Hollows
        hol_hum_depth_height[hol_hum > self.threshold_fill] = 1
        #Hummocks
        hol_hum_depth_height[hol_hum < -self.threshold_fill] = 2
        #Pure lawns
        hol_hum_depth_height[hol_hum == 0] = 0
        
        if self.number_classes is None or self.number_classes == 5:
            #lawns upper-level and lower-level
            hol_hum_depth_height[np.logical_and(hol_hum <= self.threshold_fill, hol_hum > 0)] = 3
            hol_hum_depth_height[np.logical_and(hol_hum < 0, hol_hum >= -self.threshold_fill)] = 4
            #Solving class that is filled both in DEM and inverted DEM
            hol_hum_depth_height[np.logical_and.reduce((self.hollows_lyr > 0, self.hummock_lyr > 0, self.hollows_lyr > self.hummock_lyr))] = 1
            hol_hum_depth_height[np.logical_and.reduce((self.hollows_lyr > 0, self.hummock_lyr > 0, self.hummock_lyr >= self.hollows_lyr))] = 2
            with rasterio.open(self._p_final_microtopo5, "w", **self.meta) as dest:
                dest.write(hol_hum_depth_height)
        if self.number_classes is None or self.number_classes == 3:
            #lawns upper-level and lower-level are also treated as pure lanws
            hol_hum_depth_height[np.logical_or(hol_hum_depth_height == 3, hol_hum_depth_height == 4)] = 0
            with rasterio.open(self._p_final_microtopo3, "w", **self.meta) as dest:
                dest.write(hol_hum_depth_height)
    
    def sieve_and_end(self):
        """performing sieving"""               
        if self.sieve is True and self.min_sieve_size > 0:
            print("Starting sieving...")
            if self.number_classes is None or self.number_classes == 3:
                self.perform_sieving(self._p_final_microtopo3, self._p_final_microtopo3_sieved, self.min_sieve_size)
            if self.number_classes is None or self.number_classes == 5:
                self.perform_sieving(self._p_final_microtopo5, self._p_final_microtopo5_sieved, self.min_sieve_size)
            
            final_output_path = f"{self._p_final_microtopo3_sieved} and/or {self._p_final_microtopo5_sieved}"
        else:
            final_output_path = f"{self._p_final_microtopo3} and/or {self._p_final_microtopo5}"
            
        print("The main output is saved at the path:" + final_output_path)
        print("Thanks for using the code, and validate the output if you have any ground truth")
    
    def generate_microtopography(self):
        """This method will run all methods needed to generate microtopography"""
        self.set_nodata_to_zero()
        self.fill_dem()
        self.classify_microtopography()
        self.sieve_and_end()
        
    
    @staticmethod
    def perform_sieving(input_path, output_path, min_size):
        """
        Perform sieving on a raster file to remove small features.
    
        Parameters:
            input_path (str): The path to the input raster file.
            output_path (str): The path to the output sieved raster file.
            min_size (int): The minimum number of pixels a feature must have to be retained.
    
        Returns:
            None
        """
        with rasterio.open(input_path) as src:
            # Read the raster data and metadata
            data = src.read(1, masked=True)
            profile = src.profile
    
            # Handle nodata values
            nodata_value = profile.get('nodata')
            if nodata_value is not None:
                data = np.where(data == nodata_value, -9999, data)  # Replace nodata values with a valid value
    
            # Convert data type
            data = data.astype(np.int16)
            profile['dtype'] = 'int16'
            profile['nodata'] = -9999  # Set the valid nodata value
    
            # Perform sieving
            sieved_data = sieve(data, size=min_size, connectivity=8)
    
            # Write the sieved raster data to the output file
            with rasterio.open(output_path, 'w', **profile) as dst:
                dst.write(sieved_data, 1)

    @staticmethod
    def extract_point_from_raster(path_to_raster, path_to_shapefile, colname_output):
        """Extracts raster values in a multipoint shapefile"""
        # Read points and raster files
        pts = gpd.read_file(path_to_shapefile)
        raster = rasterio.open(path_to_raster)
        #Create the column that will hold the output
        pts[colname_output] = np.nan
        #Iterating over each row of the geodataframe to extract the raster value and store it in a column in the shapefile

        for index, point in pts.iterrows():
            row, col = raster.index (point.geometry.x, point.geometry.y)
            raster_val = raster.read(1)[row, col]
            pts.loc[index, colname_output] = raster_val
            
        return pts
    
    @staticmethod
    def field_validation(path_to_raster, path_to_shapefile, colname_model, colname_field, path_out):
        """
        Performs validation given field observation as points shapefile
        
        Parameters:
            path_to_raster (str): path to the raster containing modelled values
            path_to_shapefile (str): path to the points shapefile containing ground truth
            colname_model (str): name of the column in which to store the extracted modelled values from the raster to shapefile points
            path_out (str): directory where to store the validation outputs (confusion matrix, kappa & accuracy metrics, classification report)
        """
        path_to_conf_matrix = os.path.join(path_out, "Confusion_matrix.csv")
        path_to_metrics = os.path.join(path_out, "Metrics.csv")
        new_points = HuHoLa.extract_point_from_raster(path_to_raster= path_to_raster, 
                                           path_to_shapefile = path_to_shapefile, 
                                           colname_output = colname_model)

        new_points = new_points[new_points[colname_model].notna()]
        Value = new_points[colname_model].astype(int)
        real = new_points[colname_field].astype(int)
        
        c_kappa = cohen_kappa_score(Value, real, labels=None, weights=None)
        acc = accuracy_score(Value, real)
        conf_matr = confusion_matrix(Value, real)
        classification_rep = classification_report(Value, real)

        data = [[c_kappa, acc]]
        df = pd.DataFrame(data, columns = ['Kappa', 'Accuracy'])
        
        with open(path_to_conf_matrix, 'w', newline='') as f:
            mywriter = csv.writer(f, delimiter=',')
            mywriter.writerows(conf_matr)    
            
        df.to_csv(path_to_metrics, index=False)

        with open(os.path.join(path_out, "Classification_Report.txt"), 'w') as f:
            f.write(classification_rep)

    @staticmethod
    def perform_linear_regression(df, x_col, y_col):
        """performs linear regression and calculates RMSE and r2 with equation"""
        # Extract x and y data from the DataFrame
        x = df[x_col]
        y = df[y_col]

        # Fit a linear regression model
        model = LinearRegression()
        model.fit(x.values.reshape(-1, 1), y)

        # Get the regression equation parameters
        slope = model.coef_[0]
        intercept = model.intercept_

        # Calculate R²
        y_pred = model.predict(x.values.reshape(-1, 1))
        r_squared = r2_score(y, y_pred)
        
        #calculate RMSE
        rmse = mean_squared_error(y, y_pred, squared = False)

        # Create a scatter plot with the regression line, equation, and R²
        plt.scatter(x, y, label='Data Points')
        plt.plot(x, slope * x + intercept, color='red', 
                 label=f'y = {slope:.2f}x {"+" if intercept >= 0 else "-"} {abs(intercept):.2f}\nR² = {r_squared:.2f}\nRMSE = {rmse:.2f}')
        
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.legend()
        plt.show()
    
    @staticmethod
    def convert_raster_units_m_cm(path_in, path_out):
        """Converts elevation raster units from meters to centimeters"""
        with rasterio.open(path_in) as src:
            data = src.read(1, masked=True)
            #nodata_value = src.nodata  # Get the nodata value from the original raster
            meta = src.meta
            
            data_cm = data * 100
            #meta.update(nodata=0)
            with rasterio.open(path_out, 'w', **meta) as dst:
                dst.write(data_cm, 1)
         
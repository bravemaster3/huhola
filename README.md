# HuHoLa: Hummock-Hollow-Lawn Microtopography Model
`Mires` `Peatlands` `Microtopography` `Hydromorphology` `Hydrology` `DEM` `Fill sinks`

## Description
**HuHoLa** is a straightforward model that uses digital elevation model (DEM) data to identify microtopographic features such as hummocks, hollows, and lawns. It leverages **fill algorithms** commonly used in hydrology to distinguish these features based on their topographic profiles.

- **Hollows**: Identified as depressions filled by the algorithm.
- **Hummocks**: Found by applying the fill algorithm to an inverted DEM.
- **Lawns**: Areas not filled in either the original or inverted DEM.

A more comprehensive description is available in our paper (*submission pending*).

---

## Dependencies
- **Python Libraries**: `WBT` (WhiteboxTools), `rasterio`, `numpy`
- **WhiteboxTools**: Download [here](https://www.whiteboxgeo.com/download-whiteboxtools/)

---

## Installation and Usage

### 1. Clone the Repository
Clone the repository using Git or download it manually:
- Using Git:
  ```bash
  git clone https://github.com/bravemaster3/huhola.git
  ```
- Using GitHub Desktop:
  1. Go to [HuHoLa GitHub Repository](https://github.com/bravemaster3/huhola).
  2. Click on the green "Code" button and select "Open with GitHub Desktop."
  3. Follow the instructions to clone the repository.
- Alternatively, download the repository as a zipped folder, and extract it to your preferred location.

---

### 2. Set Up the Environment
Create a Python virtual environment to manage dependencies:
- Using Conda:
  ```bash
  conda create -n huhola_env python
  conda activate huhola_env
  ```
- Alternatively, using `venv`:
  ```bash
  python -m venv huhola_env
  source huhola_env/bin/activate  # Windows: huhola_env\Scripts\activate
  ```

---

### 3. Install Dependencies
Navigate to the repository folder and install required libraries:
```bash
cd /path/to/huhola
pip install -r requirements.txt
```

---

### 4. Set Up Jupyter Notebook (Optional)
If you plan to use Jupyter notebooks:
1. Add a Jupyter kernel for the environment:
   ```bash
   python -m ipykernel install --user --name=huhola --display-name "Python (huhola)"
   ```
2. Launch JupyterLab:
   ```bash
   jupyter-lab
   ```
3. In JupyterLab, select the `Python (huhola)` kernel.

---

### 5. Download and Configure WhiteboxTools
1. Download WhiteboxTools [here](https://www.whiteboxgeo.com/download-whiteboxtools/).
2. Extract the content to a folder.
3. Note the path to the binary files; you will need it to use HuHoLa (the argument "path_wbt").

---

### 6. Explore Usage Examples
The repository includes example Jupyter notebooks in the `Examples/` folder, with a "data" subfolder containing a DEM sample to play with:
- **Microtopography Classification**: Use `huhola_microtopography.ipynb` to classify features.
- **Threshold Testing**: Try `Testing_fill_threshold.ipynb` to explore thresholds.
- **Water Table Depth Proxy**: Refer to `WTD_proxy.ipynb` for generating water table depth proxies.
  - Set `fix_flats=True` and apply a small increment to hollows.
- **Soil Temperature Proxy**: Use the hummock-hollow depth height layer from the `generate_microtopography` method.

---

## Known Issues
1. **Artificial Barriers**:
   - Obstacles may create false hummocks or dam water flow, affecting water table depth proxy accuracy.
2. **Fragmented Microtopography**:
   - Results with high-resolution DEMs (< 30 cm) can be overly fragmented. Validate your results for such data.
3. **Threshold Precision**:
   - Higher thresholds increase lawns at the expense of hummocks/hollows. Suggested starting thresholds:
     - **50 cm DEM**: 0.03–0.04
     - **30 cm DEM**: 0.05

---

## Reporting Issues
If you encounter any issues or have suggestions for improvements, please report them on the [HuHoLa GitHub Issues page](https://github.com/bravemaster3/huhola/issues). When reporting an issue:
- Describe the problem clearly.
- Include steps to reproduce the issue.
- Share any relevant error messages, screenshots, or files (if possible).

Your feedback helps make HuHoLa better!

---

## Support the Project ❤️
If you find this project helpful and would like to support further development, consider buying me a coffee:  
[![Buy Me a Coffee](https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png)](https://buymeacoffee.com/bravemaster)

Thank you for your support!

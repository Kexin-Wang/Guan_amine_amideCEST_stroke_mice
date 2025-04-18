# GuanCEST, AmineCEST, and AmideCEST Mappings of Stroke Mouse Models

This repository provides MATLAB code for generating GuanCEST, AmineCEST, and AmideCEST mappings from mouse brain MRI data. Below are the instructions for processing data acquired at 9.4T and 3T field strengths.

---

## 9.4T

To generate GuanCEST, AmineCEST, and AmideCEST mappings from the data of a **permanent** middle cerebral artery occlusion (MCAO) mouse brain:

1. **Download** all the code files from this repository.
2. **Add** the `toolbox` folder to your MATLAB path.
3. **Run** the `Main400Hz_allB1.m` script in MATLAB. This script processes the data and generates mappings based on a B₁ field strength of **1.6 μT**.

### Output
- Once the code runs, GuanCEST, amineCEST, amideCEST, and M0 images will be displayed and automatically saved in your specified folder.
- Example output image:
  ![Zguan Map Demo](9.4T/data/Zguan_map.tif)
  ![Zamine Map Demo](9.4T/data/Zamine_map.tif)
  ![Zamide Map Demo](9.4T/data/Zamide_map.tif)
  ![M0 Map Demo](9.4T/data/M0_map.tif)
---

## 3T

To generate GuanCEST and AmideCEST mappings from the data of a **transient** MCAO (tMCAO) mouse brain:

1. **Download** all the code files from this repository.
2. **Add** the `toolbox` folder to your MATLAB path.
3. **Run** the `Main3T_allB1.m` script in MATLAB. This script processes the data and generates mappings based on a B₁ field strength of **0.8 μT**.

### Output
- Once the code runs, GuanCEST, amideCEST, and M0 image will be displayed and automatically saved in your specified folder.
- Example output image:

  ![Zguan Map Demo](3T/data/Zguan_map.tif)
  ![Zamide Map Demo](3T/data/Zamide_map.tif)
  ![M0 Map Demo](3T/data/M0_map.tif)

---

## Notes
- If you encounter any issues or have questions about using the code, feel free to [open an issue](https://github.com/Kexin-Wang/Guan_amine_amideCEST_stroke_mice/issues).

---

Thank you for using this repository!




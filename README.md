# EESC 6343 – High-Resolution Target Detection and Point Cloud Map Generation for Automotive FMCW Radar

## Project Overview
This project simulates a high-resolution **Frequency-Modulated Continuous-Wave (FMCW)** radar for automotive cruise control. Using MATLAB, it models multiple moving targets, creates a virtual antenna array, and generates a 2D point cloud map.

---

## Antenna Configurations
The radar uses **TDM-MIMO** to increase the virtual aperture. We compare two setups:

### 1. High-Resolution Array (9x16)
* **Structure:** 9 Transmitters × 16 Receivers.
* **Virtual Elements:** 88 antennas.
* **Angular Resolution:** **1.3°** (Very sharp).
* **Performance:** Provides high-precision tracking with minimal errors.

### 2. Low Resolution Array (4x8)
* **Structure:** 4 Transmitters × 8 Receivers.
* **Virtual Elements:** 20 antennas.
* **Angular Resolution:** **5.7°** (Blurry).
* **Performance:** Lower precision; prone to "ghost" targets.

---

## Preprocessing Pipeline

### 1. Range-Doppler Mapping
The system processes raw radar "chirps" into a 2D grid.
* **Range:** Distance to the target.
* **Doppler:** Speed of the target relative to the radar.

### 2. Improved 2D CFAR Detection
To separate targets from noise, we use an adaptive detector:
* Calculates local noise averages.
* **Target Replacement:** Once a target is detected, it is replaced by the noise estimate in the algorithm. This ensures that strong targets do not "blind" the radar to smaller targets nearby.

### 3. Angle-of-Arrival (AoA)
By analyzing the signal across the **88 virtual elements**, we use a **Spatial FFT** to calculate the exact azimuth angle (horizontal position) of every detected car.

---

## Final Result: Point Cloud Map
The output is a spatial map comparing reality vs. radar detection:
* **X-axis:** Cross-range (Left/Right position in meters).
* **Y-axis:** Range (Distance ahead in meters).
* **Markers:** Circles represent true positions; Diamonds represent radar detections.

---

## Requirements & Execution
* **Toolboxes:** Phased Array System & Signal Processing.
* **Main Script:** `auto_phase_array.m` //used line 37 and 38 to run the simulation for both the 9x16 and 4x8 configurations (one by one).
* **Plotting:** `plot_paper_figures.m`


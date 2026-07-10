# Determinants and predictive modeling of hypoxia-induced maximal heart rate decline

This repository contains the complete codebase and analysis pipeline for the Master's Thesis in Mathematics (Università degli Studi di Trento), developed in collaboration with **CeRiSM**.

## Project Structure
* `data/`: (Optional/Placeholder) Description of the experimental datasets (Skyrunning & Multi-Study).
* `src/`: Core Python/R scripts for predictive modeling (Baseline, Regularization, Interaction models).
* `notebooks/`: Jupyter Notebooks used for Exploratory Data Analysis (EDA) and feature selection plots.
* `models/`: Saved model weights or coefficients.

## Methodology Overview
The project implements a data-driven approach to model the inter-individual variability of $HR_{max}$ reduction during acute hypoxic exposure. The modeling pipeline includes:
1. **Baseline & Regularization:** Ridge/LASSO regression to handle physiological multicollinearity.
2. **Interaction Analysis:** Multiplicative terms to capture dynamic slopes across multiple simulated altitudes.

## How to Run
1. Clone the repository: `git clone https://github.com/your-username/your-repo-name.git`
2. Install dependencies: `pip install -r requirements.txt`
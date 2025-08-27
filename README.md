
---

# 📘 Projet de Régression Bayésienne et Pénalisée en R

## Objectif

Ce projet illustre l’utilisation de différentes méthodes de **régression bayésienne et pénalisée** pour la sélection de variables et la prédiction.
L’application est réalisée sur la base de données **`telecat.csv`**, qui contient une variable réponse `Y` ainsi que des covariables explicatives (facteurs socio-culturels, préférences, etc.).

Les principales approches comparées sont :

* **RR-BLUP** (Ridge Regression Best Linear Unbiased Prediction)
* **BayesA** (modèle hiérarchique bayésien avec variance spécifique par variable)
* **LASSO bayésien (BLR)**
* **SSVS** (Stochastic Search Variable Selection)
* **Régressions pénalisées classiques (LASSO, Ridge, Elastic Net via `glmnet`)**
* **Approximate Bayesian Computation (ABC)**

---

## 📦 Installation

Avant de lancer les scripts, assurez-vous d’avoir R et RStudio installés, puis installez les packages nécessaires :

```R
install.packages(c(
  "rrBLUP", "MCMCpack", "LearnBayes", "MatrixModels", 
  "mnormt", "BLR", "glmnet"
))
```

---

##  Organisation du script

1. **Chargement des données**

   * Lecture de `telecat.csv`
   * Séparation en `X` (variables explicatives) et `Y` (variable réponse)
   * Normalisation des données
   * Split en **train (100 obs.)** et **test (50 obs.)**

2. **RR-BLUP (`rrBLUP`)**

   * Ajustement du modèle mixte
   * Estimation des effets fixes et aléatoires
   * Prédictions et corrélation préd/obs
   * Sélection des variables via les coefficients extrêmes

3. **BayesA**

   * Implémentation d’une fonction MCMC personnalisée
   * Estimation des paramètres a posteriori (`β`, `μ`, `σ²`)
   * Traces MCMC et convergence
   * Sélection de variables via les coefficients postérieurs

4. **Régression bayésienne LASSO (`BLR`)**

   * Implémentation du modèle avec pénalisation LASSO
   * Estimation et post-traitement
   * Visualisation de la distribution de `λ` (prior vs posterior)
   * Sélection des variables importantes

5. **SSVS (Stochastic Search Variable Selection)**

   * Implémentation d’une fonction SSVS maison
   * Estimation des probabilités d’inclusion des variables
   * Comparaison de plusieurs hyperparamètres (`π`, `τ0`, `τ1`)

6. **Régressions pénalisées classiques (`glmnet`)**

   * LASSO
   * Ridge
   * Elastic Net (α = 0.5)
   * Comparaison des performances en prédiction

7. **Ajout d’un effet fixe ("Sexe")**

   * Intégration dans le modèle BLR
   * Estimation de l’effet fixe et intervalle de confiance
   * Amélioration de la prédiction

8. **Approximate Bayesian Computation (ABC)**

   * Implémentation de l’algorithme ABC classique (rejet)
   * Variante avec test de Kolmogorov-Smirnov (ABC-KS)
   * Visualisation des distributions a posteriori

---

##  Résultats principaux

* **Qualité prédictive** : les différentes méthodes sont comparées via la **corrélation prédictions vs observations**.
* **Sélection de variables** : chaque méthode fournit un ensemble de variables importantes (coefficients extrêmes, inclusion SSVS, etc.).
* **Effets fixes** : l’ajout de la variable `Sexe` améliore les performances du modèle BLR.
* **ABC** : permet d’estimer des postérieurs approximatifs dans un cadre simple.

---

"# BAYESIAN-REGRESSION-R" 

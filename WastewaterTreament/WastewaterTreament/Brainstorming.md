The bioremediation process in wastewater treatment using the rotating biological contactors (RBC) is a complex and it involves both Biological and physical mechanisms.

Here, our goal is to reduce cyanide concentration while minimizng the build-up of ammonia which is a byproduct of cyanide degradayion.

Let us consider the chemical equations that governs this process:

Cyanide Degradation
$$ CN^{-}~+~H*2 O \xrightarrow{\text{Pseudomonas specie}}~ \underbrace{HCOO^{-}}*{\text{formate}}~+~NH*3$$
Ammonia Removal (Nitrification)
$$2NH_3 + O_2 \xrightarrow{\text{Nitrosomonas}}~\underbrace{2NO_2^{-}}*{\text{{intermdiate product (nitrite)}}} + ~2H^{+}~+~2H_2O$$
Nitrate to Nitrate Conversion
$$2NO_2^{-}~+~O_2\xrightarrow{\text{Nitrobacter}}~2NO_3^{-}$$

The overall chemical equation governing the remediation becomes

$$CN^{-}~+~2H_2O~+~4O_2 ~\xrightarrow{\text{Microbes}}~HCOO^{-}~+~NO_3^{-}~+~H^{+}$$

Problems:

pH Consideration - The acceptable pH for the pseudomonas species for effective degradation of cyanide falls within a neutral to slightly alkaline range i.e. 7-8, therefore the release of protons during the proposed nitrification process can lower the pH of the system and hence inhibit the performance of our main organism. We must consider the control of the pH to maintain optimal conditions for higher yield by the microbes. How do we achieve this?

## Modeling The Bioremediation Process Using Diffusion Equation

Diffusion equations are used to model the cyanide removal process in a rotating biological contactor (RBC) system, especially if spatial gradients (e.g., concentration gradients within the biofilm or reactor) are significant. Diffusion-based models are particularly useful when the system involves biofilms or heterogeneous environments where mass transfer limitations play a critical role.

### Pros & Cons

1. Th ediffusion equation captures spatial gradient and account for spatial variations in concentration (eg. within a biofilm or across the reactor), which the monod kinetics alone cannot.
2. Diffusion models explicitly include mass transfer effects, which are critical in systems where substrate (e.g., cyanide or oxygen) availability varies spatially.
3. Diffusion models are well-suited for modeling biofilms, where substrate diffusion into the biofilm and reaction kinetics are coupled.
4. In systems with non-uniform mixing or complex geometries, diffusion models provide a more accurate representation of the process.

Cons

1. Additional parameters (e.g., diffusion coefficients, biofilm thickness) are required, which may be difficult to estimate experimentally.
2. Analytical solutions are rarely possible for diffusion-reaction systems, and numerical methods are typically required.
3. Diffusion models are more complex mathematically and computationally compared to Monod kinetics-based models (computationally expensive).

- _When Do We Use The Diffusion Model_

1. When bacteria grow as biofilms on the RBC surfaces.
2. When there are significant concentration gradients within the reactor or biofilm.
3. When substrate availability is limited by diffusion rather than reaction kinetics.

### Towards The Model

- Cyanide Concentration In The Liquid Phase

$$
\dfrac{\partial C}{\partial t} = D_L \nabla^2 C - \dfrac{Q_{in}}{V}C + \dfrac{Q_{in}}{V}C_{in}-\gamma (C,X)
$$

$C$: Cyanide concentration in the liquid phase (mg/L).

$D_L$: Diffusion coefficient of cyanide in the liquid phase $(m^2/s)$

$\gamma (X,C)$ : Reaction rate of cyanide degradation $(mg/L/s)$

$\nabla^2 C$ Laplacian operator (spatial diffusion term).

- Cyanide Concentration in the biofilm:

$$
\dfrac{\partial C_b}{\partial t}= D_b\nabla^2C_b -\gamma_b(C_b,X_b)
$$

$C_b$ : Cyanide concentration in the biofilm (mg/L)

$D_b$: Diffusion coefficient of cyanide in the biofilm phase $(m^2/s)$.

$\gamma_b(C_b,X_b)$: Reaction rate of cyanide degradation in the biofilm (mg/L/s).

- Bacteria growth in the biofilm

$$
\dfrac{\partial X_b}{\partial t}=\mu(C_b)X_b -\Gamma_d X_d
$$

$X_b$: Bacteria concentration in the biofilm (mg/L).

$\mu(C_b)$: Specific groeth rate of bacteria (Monod Kinetics )
$\mu(C_b)=\mu_{max}\dfrac{C_b}{K_s + C_b}$

- Boundary Condition
- At the liquid-biofilm interface

$$
D_L \dfrac{\partial C}{\partial z}=D_b\dfrac{\partial C_b}{\partial z}
$$

- At he bottom of the biofilm (no flux):

$$
\dfrac{\partial C_b}{\partial z}=0
$$

- 1-D Model

- Liquid Phase

$$
\dfrac{\partial C}{\partial t}=D_L\dfrac{\partial^2 C}{\partial z^2}-\dfrac{Q_{in}}{V}C+\dfrac{Q_{in}}{V}C_{in}-\gamma(C,X)
$$

- Biofilm Phase

$$
\dfrac{\partial C_b}{\partial t}=D_b \dfrac{\partial^ C_b}{\partial z^2}-\Gamma_b(C_b,X_b)
$$

$$
\dfrac{\partial X_b}{\partial t}=\mu(C_b)X_b-\Gamma_d X_d
$$

- Boundary Conditions

- At $z=0$ (liquid-biofilm interface):

$$
D_L \dfrac{\partial C}{\partial t} = D_b\dfrac{\partial C_b}{\partial z}
$$

- At $z=0$

$$
\dfrac{\partial C_b}{\partial z}=0
$$

- Next step, discretize the spatial domain (e.g., using finite difference or finite element methods).

# Some Gaps To Consider

1. Limited understanding of bacteria interactions - investigate how different bacteria interact in the RBC system and how these interactions affect cyanide removal efficienc eg. How do competitive interactions impact the stability of the microbial community.
2. Can cooperative interactions (eg. cross-feeding) enhance cyanide degradation?
3. Lack of multi-dimensional modeling - develop 2 or 3D model that incoorporate spatial gradient in substrate concentration, bacteria distribution, and mass transfer. This could reveal new insight into biofilm structure and function.
4. Inadequate exploration of environmental factors. - Conduct a systematic study to quantify how variations in environmental conditions affect the performance of the bacteria species. Eg. How does temperature influence growth rate of different species and what is the optimal pH range for cyanide degradation in RBC?
5. Scaling Up from Lab to Pilot/Industrial Scale. Conduct pilot-scale studies to evaluate the feasibility of scaling up the bioremediation process. Key considerations include: How does the system perform under continuous operation? and What are the economic and logistical challenges of scaling up?
6. Long-Term Stability of Microbial Communities. Investigate the factors that influence the long-term stability of microbial communities, such as: How do operational parameters (e.g., flow rate, substrate concentration) affect community dynamics? and Can bioaugmentation (adding specific bacterial strains) help maintain community diversity?

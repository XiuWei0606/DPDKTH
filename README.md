# DPDKTH
 Simulate therma-hydraulic problems using DPDK+EDFM
 
 This simulation case DPDK_Couple_EDFM.m uses DPDK coupled with EDFM to simulate fractured reservoir.
 DPDK is used to represent tight matrix rock (m) and natural fractures (f), EDFM represents large-scale fractures.
 The grid sort order: first is f gird, then is EDFM grid, and last is m grid.

 The grid connection relationship is shown as below:
![image](https://github.com/user-attachments/assets/85cd8cc2-b04a-41c1-94d7-d5a403022fda)

 setupThermalDPDKEDFMOperatorsTPFA.m is used to get transmissibility and other information.
![image](https://github.com/user-attachments/assets/ec21fa7c-70f8-4de0-9945-6bed23bb9848)





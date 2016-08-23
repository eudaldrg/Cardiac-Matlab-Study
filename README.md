# Cardiac-Matlab-Study
Shape, motion and deformation analysis of 2D echocardiographic sequences. Examples of application to the characterization of myocardial (dys)function.

## SYNOPSIS 
 During the last years, there have been large advances on the understanding of cardiac mechanics, and in particular shape, motion and deformation. However, it is important to comprehend how these features interrelate and how they can be interpreted to refine the clinical diagnosis. The goal of this article is to illustrate the use of 2D ultrasound (US) cardiac motion sequences on some typical cardiac pathologies. Thus, a computer program has been written in MATLAB language, which allows the user to easily obtain and visualize information about myocardial shape, motion (displacement and velocity) and deformation (strain) along the left ventricle (LV) myocardium wall during the cardiac cycle. 

The method is applied to pre-processed data from 2D US sequences of four different subjects: two healthy volunteers, one patient with left ventricular mechanical dyssynchrony from left-bundle-branch-block, and one patient with hypertrophic cardiomyopathy. The interpretation based on these features was corroborated with other analyzing techniques such as electrocardiogram (ECG) and visual interpretation by experienced clinical observers.

  2D-US Sequences tracked  |  Example of the software output: Spatiotemporal representation of strain
:-------------------------:|:-------------------------:
![](https://github.com/gloriamacia/Cardiac-Matlab-Study/blob/master/Subject%20videos/02_HypertrophicCardiomyopathy-1.gif)  |  ![](https://github.com/gloriamacia/Cardiac-Matlab-Study/blob/master/html/Main_11.png)


## SETUP/USAGE/HOW TO

Download the project folder locally and open the main script (Main.m) in MATLAB. Drag the data (DATA.M) into the workspace and simply run the code. 

## SCIENTIFIC/MEDICAL/TECHNICAL INFORMATION
I provided a detailed explanation about the project in the file 
[The Heart of the Matter.pdf](https://github.com/gloriamacia/Cardiac-Matlab-Study/blob/master/The%20heart%20of%20the%20matter.pdf)

There is also a more visual summary in the file 
[summary.pdf](https://github.com/gloriamacia/Cardiac-Matlab-Study/blob/master/summary.pdf)
![alt tag](https://github.com/gloriamacia/Cardiac-Matlab-Study/blob/master/summary.png)

## LICENSE
Apache License Version 2.0 

## CONTACT / TROUBLESHOOT

This is one of my firsts contribution to open source software (still a newbie) so your feedback is really important to me. Collaborations, comments or suggestions about how to improve any aspect of the project are welcome.

If you loved the project or you have an innovative idea to ease the burden of your national health service that could be also applied to any other country... I would like to hear about it! Contact me using github or alternatively at gloriamaciamunoz@gmail.com. Before you leave, take a look at my [blog!](https://techieladyblog.wordpress.com/)

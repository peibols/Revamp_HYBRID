# arjun_HYBRID

To compile simply run
./compiler.sh main

You need to have Pythia8 installed, and make sure that the paths to the relevant folders are adapted for your computer in the compiler.sh file.

then cp main to test folder to test it.

You need to download the hydro file from this link https://drive.google.com/file/d/1ngYVF7wEuLljAC4Z6TAWs6Zdeqnx7zhR/view?usp=share_link
and put inside the test folder. This corresponds to PbPb with sqrt(s)=5.02 ATeV and centrality = 0-5%.

You also need the TAb2LL.dat (associated to initial profile), which you already have in the test folder.

Then use the script run_it.sh to run the program. There is a brief description of the meaning of the arguments in there.

Some pythia parameters are inside the setup_pythia.cmnd file, which is also already there. Note that this model is not tested with MPI=on (quenching of the different showers, hadronization...).

To use nuclear PDFs, you need to download them from here https://drive.google.com/drive/folders/1BULtpLtMBa43C-CrtSbLcuVk2mioa40L?usp=sharing. Pb is 208, and Au is 197.
We use LO for now. You will then need to put that file in the relevant pythia folder, namely .../pythia-install/share/Pythia8/xmldoc (where pythia-install is the folder where you have installed pythia).

If you do not want to use nuclear PDFs, then comment these lines out in the setup_pythia.cmnd file, by putting a "!" in front of them, such as  
!PDF:useHardNPDFB = on  
!PDF:nPDFSetB = 1  
!PDF:nPDFBeamB = 100822080  

!PDF:useHardNPDFA = on  
!PDF:nPDFSetA = 1  
!PDF:nPDFBeamA = 100822080  

The meaning of the output is the following:

...

\# event 0

Then you have info about weight, cross section, production point X and production point Y

For ex:

weight 1 cross 5.28385e-14 X 2.62739 Y 0.759099

Then the actual list of particles start. The meaning of each column is

px py pz mass pdg_id label

Each event ends with

end

...

Here are the meanings of labels for Hadron list:  

-2 = outgoing parton from hard scattering (for analysis purposes)
0 = Hadron from the jet  
1 = Positive hadron from the wake  
2 = Negative hadron from the wake

Here are the meanings of labels for Parton list:  

-2 = outgoing parton from hard scattering (for analysis purposes)  
0 = Parton from the jet

There are no wake particles at parton level.

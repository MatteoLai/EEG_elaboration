# EEG_elaboration

The files "_Patient1.mat_" and "_Patient2.mat_" contain 13 channels of single-lead EEG signals (Fs = 128 Hz) acquired on two patients during three 5-minutes phases:
- Initial relax phase (R1);
- Task: perform arithmetic operations by answering using a mouse (T);
- Final relax phase (R2);

The files "W_patient1.txt" and "W_patient2.txt" contain the demixing matrix W - respectively of patient 1 and 2 - exported after the application of ICA technique implemented by EEGLAB.

The topological maps returned by EEGLAB are saved in "_Sub1_TopologicalMap.png_" and "_Sub2_TopologicalMap.png_" - respectively of patient 1 and 2. By observing them it is possible to identify the artefacts component, and remove them to clean the signal.

In "_main_EEG_task.m_" the independent components (IC) of the signals in each channel are reconstructed, and thanks to the IC it is also possible to identify which component contributes the most to the generation of the alpha rhythm.

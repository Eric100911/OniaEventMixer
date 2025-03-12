#! /usr/bin/python3

'''
[Name of the file]
    myMixOnia.py
[Author]
    Chi Wang, Zhili College, Tsinghua University
[Description]
    Using the MixMPS module to mix the quarkonia events, producing DPS and TPS event samples from SPS samples.
    Up to this point, HELAC-Onia had provided such samples:
    - SPS: single J/psi production, single Upsilon production
    - SPS: double J/psi production, double Upsilon production, J/psi + Upsilon production
    With these samples, we can mix them to get DPS and TPS samples.
    - For instance, to obtain TPS and DPS samples for J/psi + J/psi + Upsilon,
        - DPS: SPS J/psi mixing with SPS J/psi + Upsilon; SPS J/psi + J/psi mixing with SPS Upsilon
        - TPS: SPS J/psi mixing, SPS J/psi, SPS Upsilon mixed together.
    Assigning a combination on particle level would be troubolesome, so we can use the event level to mix the events.
    - For instance, to produce DPS J/psi + Upsilon, we can create the following "recipe":
        - mixRecipe = {
        -     "sample_pp_psi_sps.lhe": 1,
        -     "sample_pp_upsilon_sps.lhe": 1
        - }
    This dictionary "mixRecipe" can be fed into the MPSEventMixer to produce the DPS J/psi + Upsilon sample.
'''

from EventUtils import EventUtils as eu
            
def main():
    mixRecipe = {
        "/eos/home-c/chiw/JpsiJpsiUps/tryHelac/HELAC-Onia-2.7.6/PROC_HO_30/P0_addon_pp_NOnia_MPS/output/sample_pp_nonia_mps.lhe": 3 # Jpsi source

        # "/eos/home-c/chiw/JpsiJpsiUps/tryHelac/HELAC-Onia-2.7.6/PROC_HO_31/P0_addon_pp_NOnia_MPS/output/sample_pp_nonia_mps.lhe": 1, # Upsilon source
    }
    eventList = eu.MPSOniaMixer.mix_from_recipe(mixRecipe, maxEvents = 50000)
    # Apply a filter to the events
    for event in eventList:
        for particle in event.particles:
            if particle.pdgId == 443 and particle.eta > 2.5:
                eventList.remove(event)
                break
    eu.LHEParser("sample_pp_psipsi_sps.lhe").setEvents(eventList).write("psipsipsi_tps_mixed.lhe")

if __name__ == "__main__":
    main()


#! /usr/bin/python3
'''
[Brief]
    This is a python script to mix three individual SPS events into one mock-up MPS event.
    This event mixing method is used to give a mock-up of TPS contribution to triple-J/psi production.
    In particular, the script is used to generate a new event with three J/psi and one gluon coming out.
[Implementation]
    The events are stored in a LHE file, which is a standard format for storing event information.
    The original files can be interpreted as gluon-gluon collision, forming a J/psi meson and a gluon.
    To ensure the validity of color flow, we pick one event as the basis and add two J/psi mesons from two other events into it.
    Solving the momentum conservation condition is not too easy. Here's a simple solution:
        - The final-state gluon will have a modified momentum to ensure the momentum conservation.
        - To preserve the gluon mass shell condition as a massless boson, we modify the gluon energy to satisfy "E = pc".
        - Up to this point, we can calculate the sum of 4-momenta of the incoming gluons.
        - In particular, the incoming gluons will not have transverse momentum.
        - This leaves us with a longitudinal momentum conservation condition and a energy conservation condition.
        - Again, using the gluon mass shell condition, we can solve the momenta of the incoming gluons.
            - Suppose we have the outgoing momenta sums up to (0, 0, pz, E).
            - The first gluon will have momentum (0, 0, (pz+E)/2, (E+pz)/2).
            - The second gluon will have momentum (0, 0, (pz-E)/2, (E-pz)/2).
[Formatting of LHE file]
    We list a few lines of the LHE file for reference:

    <LesHouchesEvents version="1.0">
    <!--
    File generated with HELAC-ONIA 
    -->
    <init>
        2212    2212  6.800000E+03  6.800000E+03       0       0   10000   10000       3       1
        9.0170635165E+07    4.1528174250E+05    1.0000000000E+00    89
    </init>
    <event>
         4    89  1.345830E+06  3.193537E+00  7.299270E-03  1.180000E-01
          21   -1    0    0  101  103  0.0000000000E+00  0.0000000000E+00  1.4763801882E-02  1.4763801882E-02  0.0000000000E+00  0.000000E+00  9.0000E+00
          21   -1    0    0  103  102  0.0000000000E+00  0.0000000000E+00 -2.8323621863E+02  2.8323621863E+02  0.0000000000E+00  0.000000E+00  9.0000E+00
         443    1    1    2    0    0 -4.8414161537E-01  6.1568455141E-01 -1.9604672526E+02  1.9607273437E+02  3.0960000000E+00  0.000000E+00  9.0000E+00
          21    1    1    2  101  102  4.8414161537E-01 -6.1568455141E-01 -8.7174729564E+01  8.7178248062E+01  0.0000000000E+00  0.000000E+00  9.0000E+00
    </event>

    In each line of particles, the info is provided as follows:
        - PDG ID, Status code, Mother1, Mother2, Color1, Color2, Px, Py, Pz, E, Mass, proper lifetime, Spin

[Usage]
    python3 myMixTPS.py sample_pp_psi_sps.lhe triple_jpsi_mixed.lhe
'''
import math

class EventUtils:

    '''
    [Name of Class]
        LorentzVector
    [Description]
        This class is used to store the 4-momentum of a particle.
    [Methods]
        __init__(self, px = 0, py = 0, pz = 0, mass = 0)
            Constructor of the class. The default value is set to 0.
        energy(self)
            Return the energy of the particle.
        momentum(self)
            Return the momentum of the particle.
        mass(self)
            Return the mass of the particle.
        set_xyzM(self, px, py, pz, mass)
            Set the 4-momentum of the particle.
        set_xyz(self, px, py, pz)
            Set the 3-momentum of the particle. The mass is kept unchanged.
        __str__(self)
            Return the string representation of the particle.
    '''
    class LorentzVector:
        P4VECTOR_TOL = 1e-4

        def __init__(self, px = 0, py = 0, pz = 0, mass = 0):
            self.__M = mass
            self.__vecXYZ = (px, py, pz)
            self.__E = (px**2 + py**2 + pz**2 + mass**2)**0.5

        @property
        def E(self):
            return self.__E

        @property
        def vecXYZ(self):
            return self.__vecXYZ

        @property
        def M(self):
            return self.__M

        @property
        def vecPerp(self):
            # As a three-vector, only with x and y components, z set to 0.
            return (self.__vecXYZ[0], self.__vecXYZ[1], 0)

        @property
        def perp(self):
            # The magnitude of the transverse component.
            return (self.__vecXYZ[0]**2 + self.__vecXYZ[1]**2)**0.5

        @property
        def eta(self):
            # Pseudorapidity of the particle.
            pz = self.__vecXYZ[2]
            return 0.5 * math.log((self.E + pz) / (self.E - pz))

        # Operator + for adding two Lorentz vectors.
        def __add__(self, other):
            res = EventUtils.LorentzVector()
            res.set_xyzE(self.__vecXYZ[0] + other.__vecXYZ[0],
                         self.__vecXYZ[1] + other.__vecXYZ[1],
                         self.__vecXYZ[2] + other.__vecXYZ[2],
                         self.__E + other.__E)
            return res

        def set_xyzM(self, px, py, pz, mass):
            self.__vecXYZ = (px, py, pz)
            self.__M = mass
            self.__E = (px**2 + py**2 + pz**2 + mass**2)**0.5
            return self

        def set_xyz(self, px, py, pz):
            self.__vecXYZ = (px, py, pz)
            self.__E = (px**2 + py**2 + pz**2 + self.__M**2)**0.5
            return self

        def set_xyzE(self, px, py, pz, energy):
            if not EventUtils.LorentzVector.is_valid(px, py, pz, 0, energy):
                raise ValueError("Mass-shell condition is not satisfied.")
            self.__vecXYZ = (px, py, pz)
            self.__E = energy
            M_square = energy**2 - (px**2 + py**2 + pz**2)
            if M_square < 0:
                self.__M = 0
            else:
                self.__M = M_square**0.5

        @staticmethod
        def is_valid(px, py, pz, mass, energy):
            return (energy**2 - (px**2 + py**2 + pz**2 + mass**2)) > ( - EventUtils.LorentzVector.P4VECTOR_TOL * (energy**2))

    '''
    [Name of Class]
        Particle
    [Description]
        This class is used to store the information of a particle.
        The info is kept consistent with the LHE file format.
    [Attributes]
        __pdgid: int
            PDG ID of the particle.
        __statuscode: int
            Status code of the particle.
        __mother1: int
            The first mother of the particle.
        __mother2: int
            The second mother of the particle.
        __color1, __color2: int
            The color flow parameters of the particle.
        __p4: LorentzVector
            The 4-momentum of the particle.
        __proper_lifetime: float
            The proper lifetime of the particle.
        __spin: float
            The cosine of the angle between the spin and the momentum.
    [Methods]
        __init__(self, pdgid, statuscode = 1, mass = 0,
                 mother1 = 0, mother2 = 0, color1 = 0, color2 = 0,
                 proper_lifetime = 0, spin = 0)
            Constructor of the class. The default value is set to 0.
        to_lhe(self)
            Return the string representation of the particle.
            Compatible with the LHE file format.
        set_momentum(self, px, py, pz)
            Set the momentum of the particle. The mass is kept unchanged.

    '''
    class Particle:
        def __init__(self, pdgid, statuscode = 1,
                     mother1 = 0, mother2 = 0, color1 = 0, color2 = 0,
                     px = 0, py = 0, pz = 0, energy = 0, mass = 0,
                     proper_lifetime = 0, spin = 0):
            # Check the mass-shell condition
            if not EventUtils.LorentzVector.is_valid(px, py, pz, mass, energy):
                raise ValueError("Mass-shell condition is not satisfied.")
            self.__pdgid = pdgid
            self.__statuscode = statuscode
            self.__mother1 = mother1
            self.__mother2 = mother2
            self.__color1 = color1
            self.__color2 = color2
            self.__p4 = EventUtils.LorentzVector(px, py, pz, mass)
            self.__proper_lifetime = proper_lifetime
            self.__spin = spin

        @property
        def pdgid(self):
            return self.__pdgid

        @property
        def statuscode(self):
            return self.__statuscode

        @property
        def mother1(self):
            return self.__mother1

        @property
        def mother2(self):
            return self.__mother2

        @property
        def color1(self):
            return self.__color1

        @property
        def color2(self):
            return self.__color2

        @property
        def p4(self):
            return self.__p4

        @property
        def proper_lifetime(self):
            return self.__proper_lifetime

        @property
        def spin(self):
            return self.__spin

        def to_lhe(self):
            return f"\t{self.__pdgid:5d}\t{self.__statuscode:3d}" + \
                   f"\t{self.__mother1:3d}\t{self.__mother2:3d}" + \
                   f"\t{self.__color1:3d}\t{self.__color2:3d}" + \
                   f"\t{self.__p4.vecXYZ[0]: 15.10E}\t{self.__p4.vecXYZ[1]: 15.10E}\t{self.__p4.vecXYZ[2]: 15.10E}\t{self.__p4.E: 15.10E}" + \
                   f"\t{self.__p4.M: 15.10E}\t{self.__proper_lifetime: 15.10E}\t{self.__spin: 15.10E}"

        def set_momentum(self, px, py, pz):
            self.__p4.set_xyz(px, py, pz)
            return self

    '''
    [Name of Class]
        Event
    [Description]
        This class is used to store the information of an event.
        The info is kept consistent with the LHE file format.
    [Attributes]
        _particles: list
            A list of particles in the event.
        __proc_id: int
            The process ID of the event.
        __scale: float
            The energy scale of the event.
        __alphaQED: float
            The QED coupling constant used to generate the event.
        __alphaQCD: float
            The QCD coupling constant used to generate the event.
        __weight: float
            The weight of the event.
    [Methods]
        __init__(self, proc_id, scale, alphaQED, alphaQCD, weight)
            Constructor of the class.
        set_from_lhe(self, lhe_block)
            Set the event information from a LHE block.
        to_lhe(self)
            Return the string representation of the event.
            Compatible with the LHE file format.
    '''
    class Event:
        STD_ALPHA_QED = 0.0072973525643
        def __init__(self, proc_id = 82, scale = 10.0, alphaQED = STD_ALPHA_QED, alphaQCD = 0.118, weight = 1.0):
            self._particles = []
            self.__proc_id = proc_id
            self.__scale = scale
            self._alphaQED = alphaQED
            self._alphaQCD = alphaQCD
            self.__weight = weight

        def set_from_lhe(self, lhe_event_block):
            # Identify the <event> and </event> tags. Remove them.
            start = lhe_event_block.index("<event>\n") + 1
            end   = lhe_event_block.index("</event>\n")
            lhe_event_block = lhe_event_block[start:end]
            # First line: number of particles, process ID, weight, scale, alphaQED, alphaQCD
            line = lhe_event_block[0].split()
            self.__proc_id    = int(line[1])
            self.__weight     = float(line[2])
            self.__scale      = float(line[3])
            self._alphaQED   = float(line[4])
            self._alphaQCD   = float(line[5])
            # Verify the number of particles
            if int(line[0]) != len(lhe_event_block) - 1:
                raise ValueError("Number of particles does not match.")
            # Parse the particle information
            self._particles = []
            for i in range(1, len(lhe_event_block)):
                # Split the line to get the particle information
                line = lhe_event_block[i].split()
                pdgid = int(line[0])
                statuscode = int(line[1])
                mother1 = int(line[2])
                mother2 = int(line[3])
                color1 = int(line[4])
                color2 = int(line[5])
                px = float(line[6])
                py = float(line[7])
                pz = float(line[8])
                energy = float(line[9])
                mass = float(line[10])
                proper_lifetime = float(line[11])
                spin = float(line[12])
                # Create a new particle object
                try:
                    particle = EventUtils.Particle(pdgid, statuscode, mother1, mother2, color1, color2,
                                                   px, py, pz, energy, mass, proper_lifetime, spin)
                    self._particles.append(particle)
                except ValueError:
                    print(f"Invalid particle: {line}")
                    continue
            return self

        def to_lhe(self):
            lhe_block = []
            lhe_block.append("<event>")
            lhe_block.append(f"\t{len(self._particles)}\t{self.__proc_id}" +
                             f"\t{self.__weight: 15.10E}\t{self.__scale: 15.10E}" +
                             f"\t{self._alphaQED: 15.10E}\t{self._alphaQCD: 15.10E}")
            for particle in self._particles:
                lhe_block.append(particle.to_lhe())
            lhe_block.append("</event>")
            return lhe_block

        @property
        def particles(self):
            return self._particles.copy()

        @property
        def nparticles(self):
            return len(self._particles)

        @property
        def proc_id(self):
            return self.__proc_id

        @property
        def scale(self):
            return self.__scale

        @property
        def alphaQED(self):
            return self._alphaQED

        @property
        def alphaQCD(self):
            return self._alphaQCD

        @property
        def weight(self):
            return self.__weight

        def count_particle(self, pdgid):
            return sum([1 for particle in self._particles if particle.pdgid == pdgid])

        def count_incoming(self):
            return sum([1 for particle in self._particles if particle.statuscode == 1])

    '''
    [Name of Class]
        OniaEvent_SPS
    [Description]
        Inherits from the Event class.
        Gluon-gluon collision event with J/psi or/and Upsilon production.
        The final state may contain at most one gluon.
    [Methods]
        __init__(self, proc_id, scale, alphaQED, alphaQCD, weight)
            Constructor of the class.

        set_from_lhe(self, lhe_block)
            Set the event information from a LHE block.
    '''
    class OniaEvent(Event):
        # Enumerate the PDG ID of the particles
        ONIA_PDGID = [
            443,        # J/psi
            100443,     # psi(2S)
            200443,     # psi(3S)
            553,        # Upsilon
            100553,     # Upsilon(2S)
            200553      # Upsilon(3S)
        ]

        # Function prototype.
        def __fix_with_gluon(self):
            pass

        # Constructor
        # The base event only provide some info. The particles are added later.
        def __init__(self, base_event, onia: list = []):
            super().__init__(base_event.proc_id, base_event.scale,
                             base_event.alphaQED, base_event.alphaQCD, base_event.weight)
            # Check if all given outgoing particles are quarkonia.
            # Also check mother particle indexing.
            self._particles = []
            for particle in onia:
                if particle.pdgid not in EventUtils.OniaEvent.ONIA_PDGID:
                    raise ValueError("Non-quarkonia particle is given.")
                if particle.mother1 != 1 or particle.mother2 != 2:
                    raise ValueError("Mother particle is not set correctly.")
                self._particles.append(particle)

            self.__fix_with_gluon()


        def __fix_with_gluon(self):
            # Add a outgoing gluon. Use its momentum to fix. Modify incoming gluons also.
            # Sum the trasnverse momentum of non-gluon outgoing particles.
            sum_p4 = EventUtils.LorentzVector()
            for particle in self._particles:
                sum_p4 += particle.p4

            # Construct a gluon to conpensate.
            new_gluon_p4 = EventUtils.LorentzVector(-sum_p4.vecXYZ[0], -sum_p4.vecXYZ[1], 0, 0)
            # Add the gluon to the list.
            new_gluon = EventUtils.Particle(pdgid=21, statuscode=1,
                                 mother1=1, mother2=2, color1=101, color2=102,
                                 px=new_gluon_p4.vecXYZ[0],
                                 py=new_gluon_p4.vecXYZ[1], 
                                 pz=new_gluon_p4.vecXYZ[2],
                                 energy=new_gluon_p4.E,
                                   mass=new_gluon_p4.M,
                                 proper_lifetime=0, spin=0)
            self._particles.append(new_gluon)

            # Modify the incoming gluons. Calculate the energy and pz.
            sum_p4 = EventUtils.LorentzVector()
            for particle in self._particles:
                sum_p4 += particle.p4

            # Colour flow settings for incoming gluons.
            # The gluon is added. Use a "triplet" colour flow.
            gluon_1_color_1 = 101
            gluon_1_color_2 = 103
            gluon_2_color_1 = 103
            gluon_2_color_2 = 102

            # Assign incoming gluon momenta.
            gluon_1_pz = (sum_p4.vecXYZ[2] + sum_p4.E) / 2
            gluon_2_pz = (sum_p4.vecXYZ[2] - sum_p4.E) / 2
            # Create the gluons.
            gluon_1 = EventUtils.Particle(pdgid=21, statuscode=-1,
                               mother1=0, mother2=0, color1=gluon_1_color_1, color2=gluon_1_color_2,
                               px=0, py=0, pz=gluon_1_pz, energy=gluon_1_pz,   mass=0,
                               proper_lifetime=0, spin=0)
            gluon_2 = EventUtils.Particle(pdgid=21, statuscode=-1,
                               mother1=0, mother2=0, color1=gluon_2_color_1, color2=gluon_2_color_2,
                               px=0, py=0, pz=gluon_2_pz, energy=-gluon_2_pz,  mass=0,
                               proper_lifetime=0, spin=0)
            self._particles.insert(0, gluon_1)
            self._particles.insert(1, gluon_2)
            return self

    '''
    [Name of class]
        MixMPS
    [Description]
        This class is used to mix SPS events into a MPS event.
    [Methods]
        mix(source_events: list)
            Mix the SPS events into a MPS event.
            The first event is used as the base event.
            The qurakonia are extracted from all events and added to form a final OniaEvent.
    '''
    class EventMixer:
        @staticmethod
        def mix(source_events: list):
            # Extract the quarkonia from all events.
            onia = []
            for event in source_events:
                for particle in event.particles:
                    if particle.pdgid in EventUtils.OniaEvent.ONIA_PDGID:
                        onia.append(particle)
            # Create a new event with the quarkonia.
            return EventUtils.OniaEvent(source_events[0], onia)

    '''
    [Name of Class]
        LHEParser
    [Description]
        This class is used to parse the LHE file. Provide input and output methods.
        An LHEParser object is hooked to a file, with read/write capability.
    [Attributes]
        __file: str
            The file name of the LHE file.
        __header: list
            The header of the LHE file, stored line by line as a list of strings.
        __Events: list
            List of events in the LHE file.
        __nEvents: int
            Number of events in the LHE file.
    [Methods]
        split_lines(file_as_list: list)
            Split the lines of the LHE file into header and event blocks.
        __init__(self, file: str)
            Constructor of the class.
            Read the LHE file and store the information.
        setEvents(self, events: list)
            Set the events of the LHE file.
        write(self, file: str)
            Write the LHE file to the given file.
        write_to_origin(self)
            Write the LHE file to the original file.
    '''
    class LHEParser:
        __file: str           # The file name of the LHE file.
        __header: list        # The header of the LHE file, stored line by line as a list of strings.
        __Events: list        # List of events in the LHE file.
        __nEvents: int        # Number of events in the LHE file.

        @staticmethod
        def split_lines(file_as_list: list):
            pass


        def __init__(self, file: str):
            self.__file = file
            self.__Events = []
            self.__nEvents = 0
            # If file already present, just read it. If not, keep this object empty.
            with open(file, 'r') as f:
                file_as_list = f.readlines()
                self.__header, events_block = EventUtils.LHEParser.split_lines(file_as_list)
                for event_block in events_block:
                    self.__Events.append(EventUtils.Event())
                    self.__Events[-1].set_from_lhe(event_block)
                self.__nEvents = len(self.__Events)

        @property
        def Events(self):
            return self.__Events.copy()

        @property
        def nEvents(self):
            return self.__nEvents

        @staticmethod
        def split_lines(file_as_list: list):
            # Use <LesHouchesEvents> and </LesHouchesEvents> as the delimiter for file. 
            # Some markers for the header and tail.
            inFile = False
            inInit = False
            inEvent = False

            # Locate the header with <LesHouchesEvents> tag.
            header_block = []
            events_block = []

            for line in file_as_list:
                if "<LesHouchesEvents" in line:
                    if not inFile:
                        inFile = True
                    else:
                        raise ValueError("Invalid LHE file format. Excess header.")

                elif not inFile:
                    continue      

                elif "</LesHouchesEvents>" in line:
                    if inFile:
                        inFile = False
                        break
                    else:
                        raise ValueError("Invalid LHE file format. Unpaired tail.") 

                elif "<init>" in line:
                    if not inInit:
                        inInit = True
                        header_block.append(line)
                    else:
                        raise ValueError("Invalid LHE file format. Excess init block.")

                elif "</init>" in line:
                    if inInit:
                        inInit = False
                        header_block.append(line)
                    else:
                        raise ValueError("Invalid LHE file format. Unpaired init end block.")

                elif inInit:
                    header_block.append(line)

                elif "<event>" in line:
                    if not inEvent:
                        inEvent = True
                        events_block.append([line])
                    else:
                        raise ValueError("Invalid LHE file format. Excess event block.")

                elif "</event>" in line:
                    if inEvent:
                        inEvent = False
                        events_block[-1].append(line)
                    else:
                        raise ValueError("Invalid LHE file format. Unpaired event end block.")

                elif inEvent:
                    events_block[-1].append(line)

            return header_block, events_block

        def setEvents(self, events: list):
            self.__Events  = events
            self.__nEvents = len(events)
            return self

        def write(self, file: str):
            with open(file, 'w') as f:
                f.write("<LesHouchesEvents version=\"1.0\">\n")
                for line in self.__header:
                    f.write(line)
                for event in self.__Events:
                    for line in event.to_lhe():
                        f.write(line + "\n")
                f.write("</LesHouchesEvents>")

        def write_to_origin(self):
            self.write(self.__file)

    class MPSOniaMixer:
        @staticmethod
        def mix_from_recipe(mixRecipe, maxEvents = 1e9):
            '''
            [Name of the function]
                mix_from_recipe
            [Description]
                Mix the quarkonia events from the files in the mixRecipe dictionary.
            [Argument]
                mixRecipe: dictionary
                    - A dictionary containing the file names and number of events mix from the file.
                    - Due to the nature of dictionary, no duplicate file names are allowed.
            [Return]
                eventList: list
                    - A list of mixed events.
            '''
            eventList = []  
            # Use eu.LHEParser to parse the LHE files that are in the mixRecipe dictionary.
            # Use a dictionary to store the parsed events, each mapped to the file name.
            fullEventDict = {}
            for fileName in mixRecipe.keys():
                fullEventDict[fileName] = EventUtils.LHEParser(fileName).Events

            # Set limit on the number of events to mix.
            for fileName in mixRecipe.keys():
                if mixRecipe[fileName] * maxEvents > len(fullEventDict[fileName]):
                    maxEvents = len(fullEventDict[fileName]) // mixRecipe[fileName]

            # Mix the events.
            for i in range(maxEvents):
                sourceEventList = []
                for fileName in mixRecipe.keys():
                    sourceEventList.extend(fullEventDict[fileName][i * mixRecipe[fileName]: (i + 1) * mixRecipe[fileName]])
                eventList.append(EventUtils.EventMixer.mix(sourceEventList)) 
            return eventList
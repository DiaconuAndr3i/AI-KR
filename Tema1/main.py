"""
Dati enter dupa fiecare solutie afisata.

Presupunem ca avem costul de mutare al unui bloc egal cu indicele in alfabet, cu indicii incepănd de la 1 (care se calculează prin 1+ diferenta dintre valoarea codului ascii al literei blocului de mutat si codul ascii al literei "a" ) . Astfel A* are trebui sa prefere drumurile in care se muta intai blocurile cu infomatie mai mica lexicografic pentru a ajunge la una dintre starile scop
"""

import copy
import os
import time
import pathlib


# informatii despre un nod din arborele de parcurgere (nu din graful initial)
class NodParcurgere:
    def __init__(self, info, parinte, nume_fisier, cost=0, h=0):
        self.info = info
        self.parinte = parinte  # parintele din arborele de parcurgere
        self.g = cost  # consider cost=1 pentru o mutare
        self.h = h
        self.f = self.g + self.h
        self.nume_fisier = nume_fisier  #Utilizam fisierul din care introducem inputul pentru a putea accesa vectorul de control self.control
        NodParcurgere.ct = 0 #Data membru utilizata pentru a obtine ordinul configuratiei in drum

    def obtineDrum(self):
        l = [self]
        nod = self
        while nod.parinte is not None:
            l.insert(0, nod.parinte)
            nod = nod.parinte
        return l

    def afisDrum(self, g, afisCost=False, afisLung=False):  # returneaza si lungimea drumului
        l = self.obtineDrum()
        for nod in l:
            #print(str(nod))
            g.write('\n')
            g.write(str(NodParcurgere.ct) + ")\n")
            NodParcurgere.ct += 1
            g.write(str(nod))
        if afisCost:
            #print("Cost: ", self.g)
            g.write("\nCost: ")
            g.write(str(self.g))
        if afisLung:
            #print("Lungime: ", len(l))
            g.write("\nLungime: ")
            g.write(str(len(l)))
        return len(l)

    def contineInDrum(self, infoNodNou):
        nodDrum = self
        while nodDrum is not None:
            if (infoNodNou == nodDrum.info):
                return True
            nodDrum = nodDrum.parinte

        return False

    def __repr__(self):
        sir = ""
        sir += str(self.info)
        return (sir)

    # euristica banală: daca nu e stare scop, returnez 1, altfel 0

    def __str__(self):
        sir = ""
        maxInalt = max([len(stiva) for stiva in self.info])
        for inalt in range(maxInalt, 0, -1):
            for stiva in self.info:
                if len(stiva) < inalt:
                    #Tinem cont de simbolurile adaugata, '(,),[,],\,/', pentru a reprezenta configuratia
                    sir += "    "
                else:
                    #Construim configuratia sub forma ceruta de problema
                    self.obj = Graph(self.nume_fisier)
                    if self.obj.control[stiva[inalt - 1]].strip() == "piramida":
                        sir += '/' + stiva[inalt - 1] + '\ '
                    elif self.obj.control[stiva[inalt - 1]].strip() == "sfera":
                        sir += '(' + stiva[inalt - 1] + ') '
                    else:
                        sir += '[' + stiva[inalt - 1] + '] '
            sir += "\n"
        # Tinem cont de simbolurile adaugata, '(,),[,],\,/', pentru a reprezenta configuratia
        sir += "#" * (4 * len(self.info) - 1)
        return sir



class Graph:  # graful problemei
    def __init__(self, nume_fisier):
        self.nume_fisier = nume_fisier
        f = open(nume_fisier, 'r')
        text_fisier = f.read()
        self.start = []
        text_fisier=text_fisier.strip().split("\n")
        self.K = text_fisier[0]     #Retinem k, numarul de stive vide
        nr_stive_configuratie = len(text_fisier[1:])        #Retinem stivele efective care apar in configuratia initiala (adica linile de la 2 pana la final)


        # Validare configuratie initiala
        for idx in range(nr_stive_configuratie):
            stiva = text_fisier[1:][idx].strip().split(",")     #Separam fieacre stiva in functie de virgula pentru a obtine 'bloculet(litera)' ca element
            # Verificam daca piramida nu se afla in varful stivei acest aspect sugerand o configuratie initiala invalida
            for item in range(len(stiva)):
                if stiva[item][0:len(stiva[item])-3] == "piramida":
                    if item+1 < len(stiva):
                        raise Exception("Configuratie initiala invalida!")
                elif stiva[item][0:len(stiva[item])-3] == "sfera":
                    # Verificam daca exista o sfera pe prima sau pe ultima stiva, existenta sugerand o configuratie invalida
                    if idx == 0:
                        if stiva[item][0:len(stiva[item])-3] == "sfera":
                            raise Exception("Configuratie initiala invalida!")
                    elif idx == nr_stive_configuratie-1:
                        if stiva[item][0:len(stiva[item])-3] == "sfera":
                            raise Exception("Configuratie initiala invalida!")
                    # Verificam sa nu existe pe pozitiile din stanga si din dreapta la acelasi nivel cu o sfera, o piramida sau spatiu liber,
                    # pentru a nu avea o configuratie invalida
                    else:
                        stivadr = text_fisier[1:][idx + 1].strip().split(",")
                        stivast = text_fisier[1:][idx - 1].strip().split(",")
                        if stiva[item][0:len(stiva[item]) - 3] == "sfera":
                            if len(stivadr) < item + 1:
                                raise Exception("Configuratie initiala invalida!")
                            if stivadr[item][0:len(stivadr[item])-3] == "piramida":
                                raise Exception("Configuratie initiala invalida!")
                            if len(stivast) < item + 1:
                                raise Exception("Configuratie initiala invalida!")
                            if stivast[item][0:len(stivast[item])-3] == "piramida":
                                raise Exception("Configuratie initiala invalida!")



        self.nr_liste = nr_stive_configuratie-int(self.K)
        configuratie = text_fisier[1:]
        self.nr_stive_vide = 0
        self.control={}
        # Construim un vector de control in functie de alfabet utilizand pentru acest fapt codurile ASCII ale literelor
        # (tinem cont de asemena si de unicitatea etichetarii bloculetelor astfel putem construi astfel contextul)
        for i in range(97,123):
            self.control[chr(i)] = ""
        for i in range(0,len(configuratie)):
            linie = configuratie[i].strip().split(",")
            per=[]
            for item in linie:
                if item == "#":
                    self.nr_stive_vide += 1
                    continue
                else:
                    per.append(item[len(item)-2])
                    self.control[item[len(item)-2]] = item[0:len(item)-3]
            # Vom obtine in self.start o lista de liste cu litere si de asemenea un vector de cntrol care va furniza natura bloculetului
            # Spre exemplu self.control['a'] va furniza 'piramida', 'cub', 'sfera' daca exista un astfel de bloculet in configuratie
            self.start.append(per)
        f.close()
        """
        self.start=[
        ['a', 'b', 'g', 'c'], 
        ['i', 'e', 'f', 'k'], 
        ['h', 'd', 'j'], 
        [], 
        ['l', 'm']
        """

    def testeaza_scop(self, nodCurent):
        # Verificam daca numarul de stive vide pentru configuratia curenta este egal cu k-ul initial,
        # in aces caz am ajuns la o configuratie finala
        stive_c = nodCurent.info
        nr_stive = len(stive_c)
        nr_stive_nevide_nod_curent=0
        for idx in range(nr_stive):
            if len(stive_c[idx]) != 0:
                nr_stive_nevide_nod_curent += 1
        if self.nr_liste == nr_stive_nevide_nod_curent:
            return 1
        return 0


        #return nodCurent.info in self.scopuri;

    # va genera succesorii sub forma de noduri in arborele de parcurgere

    def genereazaSuccesori(self, nodCurent, tip_euristica="euristica banala"):

        listaSuccesori = []
        stive_c = nodCurent.info
        nr_stive = len(stive_c)
        for idx in range(nr_stive):
            copie_interm = copy.deepcopy(stive_c)
            if len(copie_interm[idx]) == 0:
                continue
            bloc = copie_interm[idx].pop(-1)
            for j in range(nr_stive):
                stive_n = copy.deepcopy(copie_interm)
                if idx == j:
                    continue
                # Ne asiguram ca nu vom aseza niciun bloculet peste o piramida
                if len(stive_n[j]) != 0:
                    if self.control[stive_n[j][-1]].strip() == "piramida":
                        continue

                # Analizam cazul cand vrem sa mutam o sfera si omitem posibilele mutari conflictuale
                if self.control[bloc].strip() == "sfera":
                    if j == 0:
                        continue
                    elif j == nr_stive-1:
                        continue
                    elif len(stive_n[j]) + 1 > len(stive_n[j - 1]) or len(stive_n[j]) + 1 > len(stive_n[j + 1]):
                        continue
                    elif self.control[stive_n[j-1][len(stive_n[j])]].strip() == "piramida" or self.control[stive_n[j+1][len(stive_n[j])]].strip() == "piramida":
                        continue

                # Analizam cazul in care mutam un alt bloculet decat o stiva si ne asiguram ca nu generam o situatie de conflict luan un bloculet de pe o anumita pozitie
                if idx == 0:
                    if len(stive_n[idx + 1]) >= len(stive_n[idx])+1 and len(stive_n[idx + 1]) >= 1:
                        if self.control[stive_n[idx + 1][len(stive_n[idx])]].strip() == "sfera":
                            continue
                elif idx == nr_stive - 1:
                    if len(stive_n[idx - 1]) >= len(stive_n[idx])+1 and len(stive_n[idx - 1]) >= 1:
                        if self.control[stive_n[idx - 1][len(stive_n[idx])]].strip() == "sfera":
                            continue
                else:
                    if len(stive_n[idx + 1]) >= len(stive_n[idx])+1 and len(stive_n[idx + 1]) >= 1:
                        if self.control[stive_n[idx + 1][len(stive_n[idx])]].strip() == "sfera":
                            continue
                    if len(stive_n[idx - 1]) >= len(stive_n[idx])+1 and len(stive_n[idx - 1]) >= 1:
                        if self.control[stive_n[idx - 1][len(stive_n[idx])]].strip() == "sfera":
                            continue

                stive_n[j].append(bloc)

                costMutareBloc=0
                if self.control[bloc].strip() == "piramida":
                    costMutareBloc = 1
                elif self.control[bloc].strip() == "cub":
                    costMutareBloc = 2
                elif self.control[bloc].strip() == "sfera":
                    costMutareBloc = 3

                nod_nou = NodParcurgere(stive_n, nodCurent,self.nume_fisier, cost=nodCurent.g + costMutareBloc,
                                        h=self.calculeaza_h(stive_n, tip_euristica))
                if not nodCurent.contineInDrum(stive_n):
                    listaSuccesori.append(nod_nou)

        return listaSuccesori

    def calculeaza_h(self, infoNod, tip_euristica="euristica banala"):
        nr_stive = len(infoNod)
        nr_stive_nevide_nod_curent = 0
        nr_stive_vide_nod_curent = 0
        for idx in range(nr_stive):
            if len(infoNod[idx]) != 0:
                nr_stive_nevide_nod_curent += 1
            else:
                nr_stive_vide_nod_curent += 1
        # Euristica banala va returna 1 daca configuratia are un numar mai mare de stive vide decat k si 0 daca este egal
        if tip_euristica == "euristica banala":
            if self.nr_liste == nr_stive_nevide_nod_curent:
                return 0
            return 1
        # Euristica nebanala 1 numara cate stive mai trebuiesc golite pentru a obtine k stive vide si astfel pentr
        # a ajunge la o configuratie finala
        elif tip_euristica.strip() == "euristica nebanala1":
            dif = int(self.K) - nr_stive_vide_nod_curent
            return dif
        # Euristica nebanala2 porneste de la euristica nebanala 1 insa in plus cauta cele mai mici dif stive nevide pentru ca,
        # cu siguranta, numarul acesta de mutari se vor realiza spre o configuratie finala
        elif tip_euristica.strip() == "euristica nebanala2":
            dif = int(self.K) - nr_stive_vide_nod_curent
            lista = [item for item in infoNod if len(item) != 0]       # Construim o stiva auxiliara care elimina stivele vide
            lilen = list(map(lambda x: len(x), lista))       # Construim o lista de len(stiva) pentru fieacare stiva din lista
            lista_idx = sorted(range(len(lilen)), key=lambda k: lilen[k])[:dif]     # Gasim indicii celor maimici dif stive nevide
            suma = 0
            # Sumam lungimile stivelor nevide care va reprezenta numarul de mutari sigure care se realizeaza
            for item in lista_idx:
                suma += len(lista[item])
            return suma
        # Euristica neadmisibila este similara ca si mod de gandire cu euristica nebanala2 doar ca, pentru neadmisibila,
        # va aparea o supraestimare a mutarilor din presupunerea ca vom muta doar bloculete de tip sfera,
        # lucru imposibil pentru ca nu avem de exemplu in extremitati vreun astfel de bloculet
        elif tip_euristica.strip() == "euristica neadmisibila":
            dif = int(self.K) - nr_stive_vide_nod_curent
            lista = [item for item in infoNod if len(item) != 0]
            lilen = list(map(lambda x: len(x), lista))
            lista_idx = sorted(range(len(lilen)), key=lambda k: lilen[k])[:dif]
            suma = 0
            for item in lista_idx:
                suma += len(lista[item])
            return 3*suma


    def __repr__(self):
        sir = ""
        for (k, v) in self.__dict__.items():
            sir += "{} = {}\n".format(k, v)
        return (sir)



def uniform_cost(gr, timeout, g, nrSolutiiCautate, tip_euristica):
    # in coada vom avea doar noduri de tip NodParcurgere (nodurile din arborele de parcurgere)
    t1 = time.time()        # Pornim cronometrul la intrearea in functie
    c = [NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))]

    maxim = -1
    nr_noduri_calculate = 0

    while len(c) > 0:
         #calculam numarul maxim de noduri din memorie masurand lungimea cozii pentru un demers de obtinere a solutiei
        if maxim < len(c):
            maxim = len(c)

        #print("Coada actuala: " + str(c))
        #input()
        nodCurent = c.pop(0)

        if gr.testeaza_scop(nodCurent):
            t2 = time.time()
            #print("Solutie:\n", end="")
            g.write("Solutie (ucs): ")
            nodCurent.afisDrum(g, afisCost=True, afisLung=True)
            #print("\n================\n")
            g.write("\n================\n")
            #print('Timpul de gasire a unei solutii: ', round(1000 * (t2 - t1)))
            g.write('Timpul de gasire a solutiei: ')
            g.write(str(round(1000 * (t2 - t1))))       # La fieacre solutie gasita oprim inainte de afisare si astfel obtinem timpul efectiv de gasire a solutiei
            g.write(' milisec')
            g.write('\n')
            g.write('Numarul maxim de noduri existente la un moment dat in memorie: ')
            g.write(str(maxim))
            g.write('\n')
            g.write('Numar noduri calculate: ')
            g.write(str(nr_noduri_calculate))
            g.write('\n\n')
            # nr_noduri_calculate = 0
            # maxim = -1
            nrSolutiiCautate -= 1
            if nrSolutiiCautate == 0:
                return
        t2 = time.time()
        # Daca se depaseste timpul de timeout iesim dun functie si afisam un mesaj corespunzator
        if round(1000 * (t2 - t1)) > timeout:
            g.write('Timp de executie depasit pentru solutica care urmeaza!')
            return
        lSuccesori = gr.genereazaSuccesori(nodCurent, tip_euristica=tip_euristica)
        # Nodurile calculate vor fi reprezentate de fiecare succesor al nodului curent
        nr_noduri_calculate += len(lSuccesori)
        for s in lSuccesori:
            i = 0
            gasit_loc = False
            for i in range(len(c)):
                # ordonez dupa cost(notat cu g aici și în desenele de pe site)
                if c[i].g > s.g:
                    gasit_loc = True
                    break;
            if gasit_loc:
                c.insert(i, s)
            else:
                c.append(s)


# Pentru algoritmii care urmeaza, comentariile sunt similare ca la ucs



def a_star(gr, timeout, g, nrSolutiiCautate, tip_euristica):
    # in coada vom avea doar noduri de tip NodParcurgere (nodurile din arborele de parcurgere)
    t1 = time.time()
    c = [NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))]
    maxim = -1
    nr_noduri_calculate = 0
    while len(c) > 0:
        if maxim < len(c):
            maxim = len(c)
        nodCurent = c.pop(0)
        #print(nodCurent.info,'\n')
        if gr.testeaza_scop(nodCurent) == 1:
            t2 = time.time()
            g.write("Solutie (a*): ")
            #print("Solutie: ")
            nodCurent.afisDrum(g, afisCost=True, afisLung=True)
            #print("\n================\n")
            g.write("\n================\n")
            #print('Timpul de gasire a unei solutii: ', round(1000*(t2-t1)))
            g.write('Timpul de gasire a solutiei: ')
            g.write(str(round(1000 * (t2 - t1))))
            g.write(' milisec')
            g.write('\n')
            g.write('Numarul maxim de noduri existente la un moment dat in memorie: ')
            g.write(str(maxim))
            g.write('\n')
            g.write('Numar noduri calculate: ')
            g.write(str(nr_noduri_calculate))
            g.write('\n\n')
            #nr_noduri_calculate = 0
            #maxim = -1
            #input()
            nrSolutiiCautate -= 1
            if nrSolutiiCautate == 0:
                return
        t2 = time.time()
        if round(1000 * (t2 - t1)) > timeout:
            g.write('Timp de executie depasit pentru solutica care urmeaza!')
            return
        lSuccesori = gr.genereazaSuccesori(nodCurent, tip_euristica=tip_euristica)
        nr_noduri_calculate += len(lSuccesori)
        for s in lSuccesori:
            i = 0
            gasit_loc = False
            for i in range(len(c)):
                # diferenta fata de UCS e ca ordonez dupa f
                if c[i].f >= s.f:
                    gasit_loc = True
                    break;
            if gasit_loc:
                c.insert(i, s)
            else:
                c.append(s)




def a_star_optimizat(gr, timeout, g, tip_euristica):
    # in coada vom avea doar noduri de tip NodParcurgere (nodurile din arborele de parcurgere)
    t1 = time.time()
    l_open = [NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))]
    maxim = -1
    nr_noduri_calculate = 0
    # l_open contine nodurile candidate pentru expandare

    # l_closed contine nodurile expandate
    l_closed = []
    while len(l_open) > 0:
        #print("Coada actuala: " + str(l_open))
        #input()
        if maxim < len(l_open):
            maxim = len(l_open)
        nodCurent = l_open.pop(0)
        l_closed.append(nodCurent)
        if gr.testeaza_scop(nodCurent):
            t2 = time.time()
            g.write("Solutie (a*opt): ")
            #print("Solutie: \n", end="")
            nodCurent.afisDrum(g, afisCost=True, afisLung=True)
            #print("\n================\n")
            g.write("\n================\n")
            #print('Timpul de gasire a unei solutii: ', round(1000 * (t2 - t1)))
            g.write('Timpul de gasire a solutiei: ')
            g.write(str(round(1000 * (t2 - t1))))
            g.write(' milisec')
            g.write('\n')
            g.write('Numarul maxim de noduri existente la un moment dat in memorie: ')
            g.write(str(maxim))
            g.write('\n')
            g.write('Numar noduri calculate: ')
            g.write(str(nr_noduri_calculate))
            g.write('\n\n')
            # nr_noduri_calculate = 0
            # maxim = -1
            return
        t2 = time.time()
        if round(1000 * (t2 - t1)) > timeout:
            g.write('Timp de executie depasit pentru solutica care urmeaza!')
            return
        lSuccesori = gr.genereazaSuccesori(nodCurent, tip_euristica=tip_euristica)
        nr_noduri_calculate += len(lSuccesori)
        for s in lSuccesori:
            gasitC = False
            for nodC in l_open:
                if s.info == nodC.info:
                    gasitC = True
                    if s.f >= nodC.f:
                        lSuccesori.remove(s)
                    else:  # s.f<nodC.f
                        l_open.remove(nodC)
                    break

            if not gasitC:
                for nodC in l_closed:
                    if s.info == nodC.info:
                        if s.f >= nodC.f:
                            lSuccesori.remove(s)
                        else:  # s.f<nodC.f
                            l_closed.remove(nodC)
                        break
        for s in lSuccesori:
            i = 0
            gasit_loc = False
            for i in range(len(l_open)):
                # diferenta fata de UCS e ca ordonez crescator dupa f
                # daca f-urile sunt egale ordonez descrescator dupa g
                if l_open[i].f > s.f or (l_open[i].f == s.f and l_open[i].g <= s.g):
                    gasit_loc = True
                    break
            if gasit_loc:
                l_open.insert(i, s)
            else:
                l_open.append(s)


def ida_star(gr, timeout, g, nrSolutiiCautate, tip_euristica):
    t1 = time.time()
    nodStart = NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))
    nr_noduri_calculate = 0
    limita = nodStart.f
    while True:

        #print("Limita de pornire: ", limita)
        nrSolutiiCautate, rez = construieste_drum(gr, nodStart, limita, nrSolutiiCautate, tip_euristica, t1, g, nr_noduri_calculate, timeout)
        if rez == "gata":
            break
        if rez == float('inf'):
            #print("Nu exista solutii!")
            break
        limita = rez
        #print(">>> Limita noua: ", limita)
        #print('\n')


def construieste_drum(gr, nodCurent, limita, nrSolutiiCautate, tip_euristica, t1, g, nr_noduri_calculate, timeout):
    #print("A ajuns la:\n", nodCurent)
    if nodCurent.f > limita:
        return nrSolutiiCautate, nodCurent.f
    if gr.testeaza_scop(nodCurent) and nodCurent.f == limita:
        t2 = time.time()
        g.write("Solutie (ida*): ")
        #print("Solutie: ")
        nodCurent.afisDrum(g, afisCost=True, afisLung=True)
        #print("Limita la care se ajunge: ",limita)
        #print("\n================\n")
        g.write("\n================\n")
        #print('Timpul de gasire a unei solutii: ', round(1000 * (t2 - t1)))
        #input()
        g.write('Timpul de gasire a solutiei: ')
        g.write(str(round(1000 * (t2 - t1))))
        g.write(' milisec')
        g.write('\n')
        g.write('Numarul maxim de noduri existente la un moment dat in memorie: ')
        g.write(str(1))
        g.write('\n')
        g.write('Numar noduri calculate: ')
        g.write(str(nr_noduri_calculate))
        g.write('\n\n')
        nrSolutiiCautate -= 1
        if nrSolutiiCautate == 0:
            return 0, "gata"
    t2 = time.time()
    if round(1000 * (t2 - t1)) > timeout:
        g.write('Timp de executie depasit pentru solutica care urmeaza!')
        return 0, "gata"
    lSuccesori = gr.genereazaSuccesori(nodCurent, tip_euristica=tip_euristica)
    nr_noduri_calculate += len(lSuccesori)
    minim = float('inf')
    for s in lSuccesori:
        nrSolutiiCautate, rez = construieste_drum(gr, s, limita, nrSolutiiCautate, tip_euristica, t1, g, nr_noduri_calculate, timeout)
        if rez == "gata":
            return 0, "gata"
        #print("Compara ", rez, " cu ", minim)
        if rez < minim:
            minim = rez
            #print("Noul minim: ", minim)
    return nrSolutiiCautate, minim








'''
g = open('output_files\output.txt',"w+")
gr = Graph('input_files\input.txt')
timeout = 250
ida_star(gr, timeout, g, nrSolutiiCautate=15, tip_euristica="euristica banala")
g.close()
'''



print('Fisierele existente de input sunt:\n')
arr = os.listdir('./input_files')       # Afisam fisierele existente in input_files
for item in arr:
    print(item,'\n')
nume_fisier_input = input('Introduceti numele fisierului de input sub forma "input_files\example.txt":\n')
if not os.path.exists(nume_fisier_input):
    raise Exception('Fisierul introdus nu exista!')     # Aruncam o eroare daca fisierul de input nu exista

# Alegem algoritmul care urmeaza sa ruleze
print('Alegeti un algoritm din lista de mai jos:\n')
print('ucs')
print('a_star')
print('a_star_optimizat')
print('ida_star\n')
nume_algoritm = input('Introduceti numele algoritmului care va rula:\n')


# Stabilim numarul de solutii in functie de algoritm
nrsol = 0
if nume_algoritm.strip() == 'ucs' or nume_algoritm.strip() == 'a_star' or nume_algoritm.strip() == 'ida_star':
    NSOL = input('Introduceti numarul de solutii returnate de algoritm:\n')
    try:
        NSOL = int(NSOL)
    except ValueError:
        print('Nu ati introdus un numar intreg!')
    nrsol = int(NSOL)




# Stabilim timpul de timeout
timp_timeout = input('Introduceti timpul de timeout:\n')
try:
    timp_timeout = int(timp_timeout)
except ValueError:
    print('Nu ati introdus un numar intreg!')
timeout = int(timp_timeout)

# Stabilim tipul de euristica
eur = input('Introduceti tipul de euristica(ex: banala, nebanala1, nebanala2, neadmisibila):\n')
if eur.strip() != 'banala' and eur.strip() != 'nebanala1' and eur.strip() != 'nebanala2' and eur.strip() != 'neadmisibila':
    raise Exception('Tipul de euristica introdus nu exista!')
t_euristica = eur.strip()

# Generam folderul in care se vor salva output-urile daca nu exista
dirName = 'output_files'
if not os.path.exists(dirName):
    os.mkdir(dirName)
    print("Directorul", dirName,  "creat!\n")
else:
    print("Directorul", dirName,  "a fost creat in prealabil!\n")

out = os.listdir('./output_files')
print('Fisierele de output deja existente:\n')
if len(out) == 0:
    print('DIRECTOR GOL!\n')
else:
    for item in out:
        print(item, '\n')
nume_fisier_output = input('Introduceti numele fisierului de output sub forma "output_files\example.txt":\n')

if os.path.exists(nume_fisier_output.strip()):
    print('Fisierul a fost suprascris cu alta solutie!\n')

g = open(nume_fisier_output.strip(),"w+")

gr = Graph(nume_fisier_input)


# Rulam pentru fiecare algoritm in functie de parametrii introdusi
if nume_algoritm.strip() == "a_star":
    a_star(gr, timeout, g, nrSolutiiCautate=nrsol, tip_euristica='euristica ' + t_euristica)
    g.close()
    ok = 0
    # Verificam daca nu cumva starea initiala este si stare finala
    if gr.testeaza_scop(NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))):
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia initiala este si configuratie finala!')
        g.close()
        ok = 1
    # Verificam daca fisierul de output este nescris si nu este stare initiala astfel algoritmul nu a gasit o solutie si deci nu are solutii
    if os.stat(nume_fisier_output.strip()).st_size == 0 and ok == 0:
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia nu are solutii!')
        g.close()
elif nume_algoritm.strip() == "ucs":
    uniform_cost(gr, timeout, g, nrSolutiiCautate=nrsol, tip_euristica='euristica ' + t_euristica)
    g.close()
    ok = 0
    if gr.testeaza_scop(NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))):
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia initiala este si configuratie finala!')
        g.close()
        ok = 1
    if os.stat(nume_fisier_output.strip()).st_size == 0 and ok == 0:
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia nu are solutii!')
        g.close()
elif nume_algoritm.strip() == "a_star_optimizat":
    a_star_optimizat(gr,timeout, g, tip_euristica='euristica ' + t_euristica)
    g.close()
    ok = 0
    if gr.testeaza_scop(NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))):
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia initiala este si configuratie finala!')
        g.close()
        ok = 1
    if os.stat(nume_fisier_output.strip()).st_size == 0 and ok == 0:
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia nu are solutii!')
        g.close()
elif nume_algoritm.strip() == "ida_star":
    ida_star(gr, timeout, g, nrSolutiiCautate=nrsol, tip_euristica='euristica ' + t_euristica)
    g.close()
    ok = 0
    if gr.testeaza_scop(NodParcurgere(gr.start, None, gr.nume_fisier, 0, gr.calculeaza_h(gr.start))):
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia initiala este si configuratie finala!')
        g.close()
        ok = 1
    if os.stat(nume_fisier_output.strip()).st_size == 0 and ok == 0:
        g = open(nume_fisier_output.strip(), "w+")
        g.write('Configuratia nu are solutii!')
        g.close()
else:
    raise Exception('Nu ai introdus un nume de algoritm valid!')
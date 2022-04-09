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
        self.nume_fisier = nume_fisier
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
        # Selectam maximul dintre stivele unei configuratii pentru a putea construi afisarea
        maxInalt = max([len(stiva) for stiva in self.info])
        for inalt in range(maxInalt, 0, -1):
            for stiva in self.info:
                if len(stiva) < inalt:
                    sir += "  "
                else:
                    sir += stiva[inalt - 1] + " "
            sir += "\n"
            # Regula de afisare stabilita in functie de problema
        sir += "-" * (2 * len(self.info) - 1)
        return sir



class Graph:  # graful problemei
    def __init__(self, nume_fisier):
        self.nume_fisier = nume_fisier
        f = open(nume_fisier, 'r')
        text_fisier = f.read()
        self.start = []
        # Transformam intr-o lista care contine fiecare linie a fisierului de input
        text_fisier = text_fisier.strip().split("\n")
        for item in text_fisier:
            # Daca gasim vreun '=' atunci acea stiva va fi vida si trecem la iteratia urmatoare
            if item.strip() == '=':
                self.start.append([])
                continue
            # Separam dupa simbolul '#' pentru a obtine bloculet cu bloculet
            stiva = item.strip().split("#")
            stiva_parsata = []
            for elem in stiva:
                # Separam fiecare bloculet pentru a obtine id-ul sau precum si tipul sau
                tip = elem.strip().split(":")
                # Punem intr-o stiva pereche de forma (id,tip) -> ex: (4,'ps')
                stiva_parsata.append((tip[0],tip[1]))
            self.start.append(stiva_parsata)
            # Verificam ca in configuratia initiala sa nu se gaseasca vreo pisica satula
        for stiva in self.start:
            for elem in stiva:
                if elem[1].strip() == 'ps':
                    raise Exception('Starea initiala trebuie sa contina doar pisici flamande!')


        #print(self.start)

    def testeaza_scop(self, nodCurent):
        # Testam daca mai avem pisici in configuratie, daca nu atunci am gasit o stare finala
        stive = nodCurent.info
        nr_pisici_flamande = 0
        for item in stive:
            for elem in item:
                if elem[1].strip() == 'pf':
                    nr_pisici_flamande += 1
        if nr_pisici_flamande == 0:
            return 1
        return 0


        #return nodCurent.info in self.scopuri;

    # va genera succesorii sub forma de noduri in arborele de parcurgere

    def genereazaSuccesori(self, nodCurent, tip_euristica="euristica banala"):

        listaSuccesori = []
        # Accesam informatia din obiectul de tip NodParcurgere
        stiva_curenta = nodCurent.info
        nr_stive = len(stiva_curenta)
        # Parcurgem fiecare stiva din configuratia curenta
        for idx in range(nr_stive):
            # Facem o copie a configuratiei curente
            copie_initiala = copy.deepcopy(stiva_curenta)
            # Daca stiva respectiva este vida atunci nu putem lua un bloculet de pe ea si trecem la urmatoarea iteratie
            if len(copie_initiala[idx]) == 0:
                continue
            # Selectam bloculetul demutat
            bloc_de_mutat = copie_initiala[idx].pop(-1)
            # Facem o copie a configuratiei fara bloculetul pe care intentionam sa-l mutam
            copie_intermediara = copy.deepcopy(copie_initiala)
            # Parcurgem configuratia intermediara(cea fara bloculet)
            for j in range(len(copie_intermediara)):
                # Ignoram cazul cand am putea pune bloculetul din locul de unde l-am luat
                if idx == j:
                    continue
                # Facem o copie a configuratiei care va fi modificata in functie de fiecare iteratie
                copie_finala = copy.deepcopy(copie_intermediara)
                ok = 0
                # Verificam daca bloculetul pe care vrem sa-l mutam este 'pf'
                if bloc_de_mutat[1].strip() == 'pf':
                    # while-ul de mai jos este pentru bloculetele de tip 's', pentru a elimina cate un bloculet de tip
                    # 's' de pe stiva j pana cand stiva devine vida sau intalnim altceva decat 's'
                    while True:
                        # Ne asiguram ca nu incercam sa accesam un bloculet care nu exista astfel sa obtinem eroare de index
                        if len(copie_finala[j])>0:
                            if copie_finala[j][-1][1].strip() == 's':
                                # Eliminam cate un bloculet de tip 's'
                                copie_finala[j].pop()
                                # ok are rolul de a spune ca a fost facuta o mutare pentru eliminarea bloculetelor de tip 's'
                                ok = 1
                            else:
                                if ok == 1:
                                    bloc_de_mutat_copie = (bloc_de_mutat[0],'ps')
                                    copie_finala[j].append(bloc_de_mutat_copie)
                                else:
                                    copie_finala[j].append(bloc_de_mutat)
                                break

                        else:
                            copie_finala[j].append(bloc_de_mutat)
                            break
                    if ok == 0:
                        if len(copie_finala[j]) > 0:
                            # Verificam daca nu cumva punem 'pf' langa o alta 'pf' caz incompatibil
                            if copie_finala[j][-1][1].strip() == 'ps':
                                copie_finala.append(bloc_de_mutat)
                            else:
                                continue
                        else:
                            copie_finala.append(bloc_de_mutat)
                else:
                    copie_finala[j].append(bloc_de_mutat)

                # Stabilim costul de mutare
                costmutare = 0
                if bloc_de_mutat[1].strip() == 'ps':
                    costmutare = 3
                elif bloc_de_mutat[1].strip() == 'pf':
                    costmutare = 2
                elif bloc_de_mutat[1].strip() == 's':
                    costmutare = 1

                # Vom crea obiectul sa fie de tip Nodparcurgere iar ca nod si configuratie propriu-zisa va fi configuratia
                # gasita mai sus dupa mutare
                nod_nou = NodParcurgere(copie_finala, nodCurent, self.nume_fisier, cost=nodCurent.g + costmutare,
                                        h=self.calculeaza_h(copie_finala, tip_euristica))
                if not nodCurent.contineInDrum(copie_finala):
                    listaSuccesori.append(nod_nou)


        return listaSuccesori

    def calculeaza_h(self, infoNod, tip_euristica="euristica banala"):
        # Numaram cate pisisci flamande sunt in configuratie
        nr_pisici_flamande = 0
        for item in infoNod:
            for elem in item:
                if elem[1].strip() == 'pf':
                    nr_pisici_flamande += 1
        # Pentru euristica banala cu siguranta daca numarul de pisici flamande este nenul atunci
        # va ma necesita cel putin o mutare pana la starea finala
        if tip_euristica == "euristica banala":
            if nr_pisici_flamande > 0:
                return 1
            else:
                return 0
        # Pentru euristica nebanala1 nu facem altceva decat sa returnam exact numarul de pisici flamande
        # care trebuie sa devina satule
        if tip_euristica == "euristica nebanala1":
            return nr_pisici_flamande
        # Pentru euristica nebanala2 pur si simplu optimizam estimarea pentru euristica prin inmultirea
        # cu costul mutarii unei pisici flamande
        if tip_euristica == "euristica nebanala2":
            return 2 * nr_pisici_flamande
        # Pentru neadmisibila supraestimam inmultind cu cel mai mare cost posibil
        if tip_euristica == "euristica neadmisibila":
            return 3 * nr_pisici_flamande


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
        print(nodCurent.info,'\n')
        if gr.testeaza_scop(nodCurent) == 1:
            t2 = time.time()
            g.write("Solutie (a*): ")
            print("Solutie: ")
            nodCurent.afisDrum(g, afisCost=True, afisLung=True)
            print("\n================\n")
            g.write("\n================\n")
            print('Timpul de gasire a unei solutii: ', round(1000*(t2-t1)))
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







g = open('input_pisica.txt')
gr = Graph('input_pisica.txt')
timeout = 10000
nrSolutiiCautate = 1
a_star(gr, timeout, g, nrSolutiiCautate, tip_euristica="euristica banala")
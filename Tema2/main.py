import time
import copy
import pygame
import sys
from statistics import median

ADANCIME_MAX = 1
VECT_NR_TOTAL = []

def elem_identice(lista):
    if (len(set(lista)) == 1):
        return lista[0] if lista[0] != Joc.GOL else False
    return False


class Joc:
    """
    Clasa care defineste jocul. Se va schimba de la un joc la altul.
    """
    JMIN = None
    JMAX = None
    GOL = '#'
    NR_LINII = None
    NR_COLOANE = None
    scor_maxim = 0

    # Pentru rerinerea scorului celor 2 oponenti
    scor_JMIN = 0
    scor_JMAX = 0

    def __init__(self, matr=None, NR_LINII=None, NR_COLOANE=None):
        # creez proprietatea ultima_mutare # (l,c)
        self.ultima_mutare = None
        #Adaugam penultima mutare pentru a verifica daca aduce sau nu punct
        self.penultima_mutare = None

        if matr:
            # e data tabla, deci suntem in timpul jocului
            self.matr = matr
        else:
            # nu e data tabla deci suntem la initializare
            self.matr = [[self.__class__.GOL] * NR_COLOANE for i in range(NR_LINII)]
            self.matr[4][4] = 'x'
            self.matr[4][5] = '0'
            self.matr[5][4] = 'x'
            self.matr[5][5] = '0'

            if NR_LINII is not None:
                self.__class__.NR_LINII = NR_LINII
            if NR_COLOANE is not None:
                self.__class__.NR_COLOANE = NR_COLOANE

            ######## calculare scor maxim ###########
            #sc_randuri = (NR_COLOANE - 3) * NR_LINII
            #sc_coloane = (NR_LINII - 3) * NR_COLOANE
            sc_diagonale = (NR_LINII - 3) * (NR_COLOANE - 3) * 2
            self.__class__.scor_maxim = sc_diagonale
            #Functia deseneaza_castigator coloreaza celula de pe linia i coloana j
    def deseneaza_castigator(self, lin=None, col=None):
        culoare = (0, 191, 255)
        for ind in range(self.__class__.NR_COLOANE * self.__class__.NR_LINII):
            linie = ind // self.__class__.NR_COLOANE  # // inseamna div
            coloana = ind % self.__class__.NR_COLOANE
            if linie == lin and coloana == col:
                pygame.draw.rect(self.__class__.display, culoare, self.__class__.celuleGrid[ind])
                if self.matr[linie][coloana] == 'x':
                    self.__class__.display.blit(self.__class__.x_img, (
                        coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
                elif self.matr[linie][coloana] == '0':
                    self.__class__.display.blit(self.__class__.zero_img, (
                        coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
            # pygame.display.flip()
        pygame.display.update()

    def deseneaza_grid(self, coloana_marcaj=None):  # tabla de exemplu este ["#","x","#","0",......]

        for ind in range(self.__class__.NR_COLOANE * self.__class__.NR_LINII):
            linie = ind // self.__class__.NR_COLOANE  # // inseamna div
            coloana = ind % self.__class__.NR_COLOANE

            if coloana == coloana_marcaj:
                # daca am o patratica selectata, o desenez cu rosu
                culoare = (255, 255, 0)
            else:
                # altfel o desenez cu alb
                culoare = (255, 255, 255)
            pygame.draw.rect(self.__class__.display, culoare, self.__class__.celuleGrid[ind])  # alb = (255,255,255)
            if self.matr[linie][coloana] == 'x':
                self.__class__.display.blit(self.__class__.x_img, (
                    coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
            elif self.matr[linie][coloana] == '0':
                self.__class__.display.blit(self.__class__.zero_img, (
                    coloana * (self.__class__.dim_celula + 1), linie * (self.__class__.dim_celula + 1)))
        # pygame.display.flip()
        pygame.display.update()

    @classmethod
    def jucator_opus(cls, jucator):
        return cls.JMAX if jucator == cls.JMIN else cls.JMIN

    @classmethod
    def initializeaza(cls, display, NR_LINII=6, NR_COLOANE=7, dim_celula=100):
        cls.display = display
        cls.dim_celula = dim_celula
        cls.x_img = pygame.image.load('ics.png')
        cls.x_img = pygame.transform.scale(cls.x_img, (dim_celula, dim_celula))
        cls.zero_img = pygame.image.load('zero.png')
        cls.zero_img = pygame.transform.scale(cls.zero_img, (dim_celula, dim_celula))
        cls.celuleGrid = []  # este lista cu patratelele din grid
        for linie in range(NR_LINII):
            for coloana in range(NR_COLOANE):
                patr = pygame.Rect(coloana * (dim_celula + 1), linie * (dim_celula + 1), dim_celula, dim_celula)
                cls.celuleGrid.append(patr)

    def parcurgere(self, directie):
        um = self.ultima_mutare  # (l,c)
        culoare = self.matr[um[0]][um[1]]
        nr_mutari = 0
        while True:
            um = (um[0] + directie[0], um[1] + directie[1])
            if not 0 <= um[0] < self.__class__.NR_LINII or not 0 <= um[1] < self.__class__.NR_COLOANE:
                break
            if not self.matr[um[0]][um[1]] == culoare:
                break
            nr_mutari += 1
        return nr_mutari

    # Functia anteparcurgere este similiara functiei parcurgere doar ca,
    # se va trata cazul penultimei mutari pentru a verifica daca se formeaza o diagonala
    # de 3 simboluri ale jucatorului curent
    def anteparcurgere(self, directie):
        pm = self.penultima_mutare  # (l,c)
        culoare = self.matr[pm[0]][pm[1]]
        nr_mutari = 0
        while True:
            pm = (pm[0] + directie[0], pm[1] + directie[1])
            if not 0 <= pm[0] < self.__class__.NR_LINII or not 0 <= pm[1] < self.__class__.NR_COLOANE:
                break
            if not self.matr[pm[0]][pm[1]] == culoare:
                break
            nr_mutari += 1
        return nr_mutari

    def final(self):
        if not self.ultima_mutare:
            return False
        # directii = [[(0, 1), (0, -1)], [(1, 1), (-1, -1)], [(1, -1), (-1, 1)], [(1, 0), (-1, 0)]]
        # Vectorul de directii pentru diagonale
        directii = [[(1, 1), (-1, -1)], [(1, -1), (-1, 1)]]
        um = self.ultima_mutare
        rez = False
        for per_dir in directii:
            len_culoare = self.parcurgere(per_dir[0]) + self.parcurgere(per_dir[1]) + 1  # +1 pt chiar ultima mutare
            if len_culoare >= 3:
                rez = self.matr[um[0]][um[1]]
        if (rez):
            return rez
        elif all(self.__class__.GOL not in x for x in self.matr):
            return 'remiza'
        else:
            return False

    # functia antefinal este similara functiei final insa aceasta lucreaza
    # in functie de penultima mutare, aspect care ne permite verificarea obtinerii
    # unui punt pentru penultima mutare
    def antefinal(self):
        if not self.penultima_mutare:
            return False
        # directii = [[(0, 1), (0, -1)], [(1, 1), (-1, -1)], [(1, -1), (-1, 1)], [(1, 0), (-1, 0)]]
        # Vectorul de directii pentru diagonale
        directii = [[(1, 1), (-1, -1)], [(1, -1), (-1, 1)]]
        pm = self.penultima_mutare
        rez = False
        for per_dir in directii:
            len_culoare = self.anteparcurgere(per_dir[0]) + self.anteparcurgere(per_dir[1]) + 1  # +1 pt chiar ultima mutare
            if len_culoare >= 3:
                rez = self.matr[pm[0]][pm[1]]
        if (rez):
            return rez
        elif all(self.__class__.GOL not in x for x in self.matr):
            return 'remiza'
        else:
            return False


    def mutari(self, jucator):
        l_mutari = []
        for i in range(self.__class__.NR_LINII):
            for j in range(self.__class__.NR_COLOANE):
                # Dca celula nu este libera trecem la urmatoarea
                if self.matr[i][j] != self.__class__.GOL:
                    continue
                # Facem o copie si plasam prima piesa din placuta
                matr_tabla_noua = copy.deepcopy(self.matr)
                matr_tabla_noua[i][j] = jucator
                # Stabilim un vector de directii corespunzator verticalei si orizontalei,
                # care ne va permite plasarea celei de-a doua piese astfel incat sa
                # se formeze o placuta de 1x2 sau 2x1
                dir = [[0, -1], [-1, 0], [0, 1], [1, 0]]
                for item in dir:
                    # Verificam sa nu depasim extremitatile grid-ului
                    if i + item[0] < 0 or j + item[1] < 0 or i + item[0] > self.__class__.NR_LINII - 1 or j + \
                            item[1] > self.__class__.NR_COLOANE - 1:
                        continue
                    # Fcaem o copie si plasam a doua piesa din placuta daca se poate
                    copie = copy.deepcopy(matr_tabla_noua)
                    if copie[i + item[0]][j + item[1]] == self.__class__.GOL:
                        copie[i + item[0]][j + item[1]] = jucator
                    else:
                        continue


                    # Ne asiguram ca am construit placuta pe aceasi linie pentru a avea un
                    # reper in verificarea vecinilor
                    if item[0] == 0: #Pe aceeasi linie
                        nrx = 0
                        nr0 = 0
                        # Luam minimul dupa coloana pentru a stabili unde ne situam
                        # in raport cu placuta de 1x2 plasata
                        minim = min(j, j + item[1])
                        # Consideram un vector de directii care contine toti vecinii unei placute de 1x2
                        # pentru care consideram vecinii in raport cu minimul coloanei
                        dir_ver = [[0, -1], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [0, 2], [1, 2], [1, 1], [1, 0], [1, -1]]
                        for elem in dir_ver:
                            # Verificam sa nu depasim limitele gridului
                            if i + elem[0] < 0 or minim + elem[1] < 0 or i + elem[0] > self.__class__.NR_LINII - 1 or minim + \
                                    elem[1] > self.__class__.NR_COLOANE - 1:
                                continue
                            # Contorizam numarul de x si 0 din jurul placutei
                            if copie[i + elem[0]][minim + elem[1]] == 'x':
                                nrx += 1
                            elif copie[i + elem[0]][minim + elem[1]] == '0':
                                nr0 += 1
                        # Daca avem cel putin un x si un 0 atunci este o mutare valida si adaugam
                        # in lisat de mutari
                        if nrx >= 1 and nr0 >= 1:
                            jn = Joc(copie)
                            lista = [item.matr for item in l_mutari]
                            if jn.matr not in lista:
                                jn.ultima_mutare = (i + item[0], j + item[1])
                                jn.penultima_mutare = (i,j)
                                l_mutari.append(jn)

                    # Procedam asemenator si pentru cazul in care ne situam pe aceasi coloana
                    # si consideram de aceasta data minimul pe linii si deci avem o placuta de 2x1
                    elif item[1] == 0: #Pe aceeasi coloana
                        nrx = 0
                        nr0 = 0
                        minim = min(i, i + item[0])
                        dir_ver = [[0, -1], [-1, -1], [-1, 0], [-1, 1], [0, 1], [1, 1], [2, 1], [2, 0], [2, -1],[1, -1]]
                        for elem in dir_ver:
                            if minim + elem[0] < 0 or j + elem[1] < 0 or minim + elem[
                                0] > self.__class__.NR_LINII - 1 or j + \
                                    elem[1] > self.__class__.NR_COLOANE - 1:
                                continue
                            if copie[minim + elem[0]][j + elem[1]] == 'x':
                                nrx += 1
                            elif copie[minim + elem[0]][j + elem[1]] == '0':
                                nr0 += 1
                        if nrx >= 1 and nr0 >= 1:
                            jn = Joc(copie)
                            lista = [item.matr for item in l_mutari]
                            if jn.matr not in lista:
                                jn.ultima_mutare = (i + item[0], j + item[1])
                                jn.penultima_mutare = (i,j)
                                l_mutari.append(jn)

        lista = []


        #############################################
        # Sortarea

        # Stabilim obiectele sa fie de tip Stare pentru a putea accesa scorul
        # si functia de estimare si construim o lista de tupluri pentru a putea
        # selecta usor obiectele initiale de tip Joc
        for item in l_mutari:
            stare_curenta = Stare(item, jucator, adancime=1)
            stare_curenta.scor = stare_curenta.tabla_joc.estimeaza_scor(stare_curenta.adancime)
            lista.append([stare_curenta, stare_curenta.scor, item])
        # Sortam lista de tupluri dupa al doilea parametru, dupa scor
        l_sort = sorted(lista, key=lambda x: int(x[1]))
        l_sort_final = []

        # Extragem din lista de tupluri doar acele obiecte de tip Joc
        # ordonate de aceasta data dupa scor
        for i in range(0,len(l_sort)):
            l_sort_final.append(l_sort[i][2])

        return l_sort_final


    # linie deschisa inseamna linie pe care jucatorul mai poate forma o configuratie castigatoare
    # practic e o linie fara simboluri ale jucatorului opus
    def linie_deschisa(self, lista, jucator, estimare_scor_tip = 'tip2'):
        jo = self.jucator_opus(jucator)
        # verific daca pe linia data nu am simbolul jucatorului opus
        if not jo in lista:
            if estimare_scor_tip == 'tip1':
                return 1

            # Voi returna numarul de simboluri ale jucatorului curent
            if estimare_scor_tip == 'tip2':
                return lista.count(jucator)
        return 0

    def linii_deschise(self, jucator, estimare_scor_tip='tip2'):

        linii = 0


        # Am pastrat doar cazurile de linii deschise pentru diagonale
        # si am adaptat algoritmul la 3 simboluri succesive in diagonala
        # \
        for i in range(self.__class__.NR_LINII - 2):
            for j in range(self.__class__.NR_COLOANE - 2):
                linii += self.linie_deschisa([self.matr[i + k][j + k] for k in range(0, 3)], jucator, estimare_scor_tip=estimare_scor_tip)

        # /
        for i in range(self.__class__.NR_LINII - 2):
            for j in range(2, self.__class__.NR_COLOANE):
                linii += self.linie_deschisa([self.matr[i + k][j - k] for k in range(0, 3)], jucator, estimare_scor_tip=estimare_scor_tip)

        return linii

    def estimeaza_scor(self, adancime, estimare_scor_tip='tip2'):
        t_final = self.final()
        # if (adancime==0):
        if t_final == self.__class__.JMAX:
            return (self.__class__.scor_maxim + adancime)
        elif t_final == self.__class__.JMIN:
            return (-self.__class__.scor_maxim - adancime)
        elif t_final == 'remiza':
            return 0
        else:
            return (self.linii_deschise(self.__class__.JMAX,estimare_scor_tip=estimare_scor_tip) - self.linii_deschise(self.__class__.JMIN,estimare_scor_tip=estimare_scor_tip))

    def sirAfisare(self):
        sir = "  |"
        sir += " ".join([str(i) for i in range(self.NR_COLOANE)]) + "\n"
        sir += "-" * (self.NR_COLOANE + 1) * 2 + "\n"
        sir += "\n".join([str(i) + " |" + " ".join([str(x) for x in self.matr[i]]) for i in range(len(self.matr))])
        return sir

    def __str__(self):
        return self.sirAfisare()

    def __repr__(self):
        return self.sirAfisare()


class Stare:
    """
    Clasa folosita de algoritmii minimax si alpha-beta
    Are ca proprietate tabla de joc
    Functioneaza cu conditia ca in cadrul clasei Joc sa fie definiti JMIN si JMAX (cei doi jucatori posibili)
    De asemenea cere ca in clasa Joc sa fie definita si o metoda numita mutari() care ofera lista cu configuratiile posibile in urma mutarii unui jucator
    """

    def __init__(self, tabla_joc, j_curent, adancime, parinte=None, scor=None):
        self.tabla_joc = tabla_joc
        self.j_curent = j_curent

        # adancimea in arborele de stari
        self.adancime = adancime

        # scorul starii (daca e finala) sau al celei mai bune stari-fiice (pentru jucatorul curent)
        self.scor = scor

        # lista de mutari posibile din starea curenta
        self.mutari_posibile = []

        # cea mai buna mutare din lista de mutari posibile pentru jucatorul curent
        self.stare_aleasa = None

    def mutari(self):
        l_mutari = self.tabla_joc.mutari(self.j_curent)
        juc_opus = Joc.jucator_opus(self.j_curent)
        l_stari_mutari = [Stare(mutare, juc_opus, self.adancime - 1, parinte=self) for mutare in l_mutari]

        return l_stari_mutari

    def __str__(self):
        sir = str(self.tabla_joc) + "\n(Juc curent:" + self.j_curent + ")\n"
        return sir

    def __repr__(self):
        sir = str(self.tabla_joc) + "(Juc curent:" + self.j_curent + ")\n"
        return sir


""" Algoritmul MinMax """


def min_max(stare,estimare_scor_tip='tip2'):

    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime,estimare_scor_tip=estimare_scor_tip)
        return stare

    # calculez toate mutarile posibile din starea curenta
    global nr_noduri_generate
    nr_noduri_generate = 0
    stare.mutari_posibile = stare.mutari()
    nr_noduri_generate = len(stare.mutari_posibile)
    VECT_NR_TOTAL.append(nr_noduri_generate)


    # aplic algoritmul minimax pe toate mutarile posibile (calculand astfel subarborii lor)

    mutari_scor = [min_max(mutare,estimare_scor_tip=estimare_scor_tip) for mutare in stare.mutari_posibile]

    if stare.j_curent == Joc.JMAX:
        # daca jucatorul e JMAX aleg starea-fiica cu scorul maxim
        stare.stare_aleasa = max(mutari_scor, key=lambda x: x.scor)
    else:
        # daca jucatorul e JMIN aleg starea-fiica cu scorul minim
        stare.stare_aleasa = min(mutari_scor, key=lambda x: x.scor)
    stare.scor = stare.stare_aleasa.scor
    return stare


# Functia min_max_generator este utila in momentele in care daca placuta de 1x2 sau 2x1
# aduce 2 puncte, va genera in continuare mutarea calculatorului si nu va bloca la
# primul return din functia min_max de mai sus
def min_max_generator(stare,estimare_scor_tip='tip2'):

    # calculez toate mutarile posibile din starea curenta
    stare.mutari_posibile = stare.mutari()

    # aplic algoritmul minimax pe toate mutarile posibile (calculand astfel subarborii lor)
    mutari_scor = [min_max(mutare,estimare_scor_tip=estimare_scor_tip) for mutare in stare.mutari_posibile]

    if stare.j_curent == Joc.JMAX:
        # daca jucatorul e JMAX aleg starea-fiica cu scorul maxim
        stare.stare_aleasa = max(mutari_scor, key=lambda x: x.scor)
    else:
        # daca jucatorul e JMIN aleg starea-fiica cu scorul minim
        stare.stare_aleasa = min(mutari_scor, key=lambda x: x.scor)
    stare.scor = stare.stare_aleasa.scor
    return stare


def alpha_beta(alpha, beta, stare,estimare_scor_tip='tip2'):
    if stare.adancime == 0 or stare.tabla_joc.final():
        stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime,estimare_scor_tip=estimare_scor_tip)
        return stare
    if alpha > beta:
        return stare  # este intr-un interval invalid deci nu o mai procesez

    global nr_noduri_generate
    #nr_noduri_generate = 0
    stare.mutari_posibile = stare.mutari()
    nr_noduri_generate = len(stare.mutari_posibile)
    VECT_NR_TOTAL.append(nr_noduri_generate)

    if stare.j_curent == Joc.JMAX:
        scor_curent = float('-inf')

        for mutare in stare.mutari_posibile:
            # calculeaza scorul
            stare_noua = alpha_beta(alpha, beta, mutare,estimare_scor_tip=estimare_scor_tip)

            if (scor_curent < stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor
            if (alpha < stare_noua.scor):
                alpha = stare_noua.scor
                if alpha >= beta:
                    break

    elif stare.j_curent == Joc.JMIN:
        scor_curent = float('inf')

        for mutare in stare.mutari_posibile:

            stare_noua = alpha_beta(alpha, beta, mutare,estimare_scor_tip=estimare_scor_tip)

            if (scor_curent > stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor

            if (beta > stare_noua.scor):
                beta = stare_noua.scor
                if alpha >= beta:
                    break
    stare.scor = stare.stare_aleasa.scor
    return stare

# Este similara constructiei lui min_max_generator
def alpha_beta_generator(alpha, beta, stare,estimare_scor_tip='tip2'):
    if alpha > beta:
        return stare  # este intr-un interval invalid deci nu o mai procesez

    stare.mutari_posibile = stare.mutari()


    if stare.j_curent == Joc.JMAX:
        scor_curent = float('-inf')

        for mutare in stare.mutari_posibile:
            # calculeaza scorul
            stare_noua = alpha_beta(alpha, beta, mutare,estimare_scor_tip=estimare_scor_tip)

            if (scor_curent < stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor
            if (alpha < stare_noua.scor):
                alpha = stare_noua.scor
                if alpha >= beta:
                    break

    elif stare.j_curent == Joc.JMIN:
        scor_curent = float('inf')

        for mutare in stare.mutari_posibile:

            stare_noua = alpha_beta(alpha, beta, mutare,estimare_scor_tip=estimare_scor_tip)

            if (scor_curent > stare_noua.scor):
                stare.stare_aleasa = stare_noua
                scor_curent = stare_noua.scor

            if (beta > stare_noua.scor):
                beta = stare_noua.scor
                if alpha >= beta:
                    break
    stare.scor = stare.stare_aleasa.scor
    return stare


def afis_daca_final(stare_curenta):
    final = stare_curenta.tabla_joc.final()
    if (final):
        return True
    return False

# Functia este similara lui afis_daca_final doar ca va lucra pe
# functiile mentionate mai sus, antefinal si anteparcurgere in acest fel prelucrand si penultima mutare
def afis_daca_final_pc(stare_curenta):
    antefinal = stare_curenta.tabla_joc.antefinal()
    if (antefinal):
        return True
    return False

class Buton:
    def __init__(self, display=None, left=0, top=0, w=0, h=0, culoareFundal=(53, 80, 115),
                 culoareFundalSel=(89, 134, 194), text="", font="arial", fontDimensiune=16, culoareText=(255, 255, 255),
                 valoare=""):
        self.display = display
        self.culoareFundal = culoareFundal
        self.culoareFundalSel = culoareFundalSel
        self.text = text
        self.font = font
        self.w = w
        self.h = h
        self.selectat = False
        self.fontDimensiune = fontDimensiune
        self.culoareText = culoareText
        # creez obiectul font
        fontObj = pygame.font.SysFont(self.font, self.fontDimensiune)
        self.textRandat = fontObj.render(self.text, True, self.culoareText)
        self.dreptunghi = pygame.Rect(left, top, w, h)
        # aici centram textul
        self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)
        self.valoare = valoare

    def selecteaza(self, sel):
        self.selectat = sel
        self.deseneaza()

    def selecteazaDupacoord(self, coord):
        if self.dreptunghi.collidepoint(coord):
            self.selecteaza(True)
            return True
        return False

    def updateDreptunghi(self):
        self.dreptunghi.left = self.left
        self.dreptunghi.top = self.top
        self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)

    def deseneaza(self):
        culoareF = self.culoareFundalSel if self.selectat else self.culoareFundal
        pygame.draw.rect(self.display, culoareF, self.dreptunghi)
        self.display.blit(self.textRandat, self.dreptunghiText)


class GrupButoane:
    def __init__(self, listaButoane=[], indiceSelectat=0, spatiuButoane=10, left=0, top=0):
        self.listaButoane = listaButoane
        self.indiceSelectat = indiceSelectat
        self.listaButoane[self.indiceSelectat].selectat = True
        self.top = top
        self.left = left
        leftCurent = self.left
        for b in self.listaButoane:
            b.top = self.top
            b.left = leftCurent
            b.updateDreptunghi()
            leftCurent += (spatiuButoane + b.w)

    def selecteazaDupacoord(self, coord):
        for ib, b in enumerate(self.listaButoane):
            if b.selecteazaDupacoord(coord):
                self.listaButoane[self.indiceSelectat].selecteaza(False)
                self.indiceSelectat = ib
                return True
        return False

    def deseneaza(self):
        # atentie, nu face wrap
        for b in self.listaButoane:
            b.deseneaza()

    def getValoare(self):
        return self.listaButoane[self.indiceSelectat].valoare


############# ecran initial ########################



################################################
# Functia va permite alegerea simbolului x sau 0
def x_0(display):
    btn_x_0 = GrupButoane(
        top=120,
        left=165,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="x", valoare="x"),
            Buton(display=display, w=80, h=30, text="0", valoare="0"),
        ],
        indiceSelectat=0)
    ok = Buton(display=display, top=200, left=230, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_x_0.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_x_0.selecteazaDupacoord(pos):
                    if ok.selecteazaDupacoord(pos):
                        display.fill((0, 0, 0))  # stergere ecran
                        return btn_x_0.getValoare()
        pygame.display.update()



################################################
# Functia va permite alegerea algoritmului de rulare
def algoritm(display):
    btn_algoritm = GrupButoane(
        top=120,
        left=165,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="minimax", valoare="minimax"),
            Buton(display=display, w=80, h=30, text="alphabeta", valoare="alphabeta"),
        ],
        indiceSelectat=1)
    ok = Buton(display=display, top=200, left=230, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_algoritm.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_algoritm.selecteazaDupacoord(pos):
                    if ok.selecteazaDupacoord(pos):
                        display.fill((0, 0, 0))  # stergere ecran
                        return btn_algoritm.getValoare()
        pygame.display.update()



###############################################
# Functia va permite alegerea tipului de estimare a scorului
def tip_estimare(display):
    btn_estimare = GrupButoane(
        top=120,
        left=165,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="tip1", valoare="tip1"),
            Buton(display=display, w=80, h=30, text="tip2", valoare="tip2"),
        ],
        indiceSelectat=1)
    ok = Buton(display=display, top=200, left=230, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_estimare.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_estimare.selecteazaDupacoord(pos):
                    if ok.selecteazaDupacoord(pos):
                        display.fill((0, 0, 0))  # stergere ecran
                        return btn_estimare.getValoare()
        pygame.display.update()


##########################################
# Alegerea dificultatii
def dificultate(display):
    btn_dificultate = GrupButoane(
        top=120,
        left=110,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="scazuta", valoare="scazuta"),
            Buton(display=display, w=80, h=30, text="medie", valoare="medie"),
            Buton(display=display, w=100, h=30, text="ridicata", valoare="ridicata")
        ],
        indiceSelectat=0)
    ok = Buton(display=display, top=200, left=230, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_dificultate.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_dificultate.selecteazaDupacoord(pos):
                    if ok.selecteazaDupacoord(pos):
                        display.fill((0, 0, 0))  # stergere ecran
                        return btn_dificultate.getValoare()
        pygame.display.update()



####################################################
# Alegerea tipului de joc
def oponenti(display):
    btn_oponenti = GrupButoane(
        top=120,
        left=110,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="1v1", valoare="1v1"),
            Buton(display=display, w=80, h=30, text="vs PC", valoare="cmpvspl"),
            Buton(display=display, w=100, h=30, text="PC1 vs PC2", valoare="pc1vspc2")
        ],
        indiceSelectat=1)
    ok = Buton(display=display, top=200, left=230, w=40, h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_oponenti.deseneaza()
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not btn_oponenti.selecteazaDupacoord(pos):
                    if ok.selecteazaDupacoord(pos):
                        display.fill((0, 0, 0))  # stergere ecran
                        return btn_oponenti.getValoare()
        pygame.display.update()



#############################################################
# Functia de afisare principala care va apela succesiv functiile anterior
# mentionate si care va permite jocului sa porneasca
def incepe_joc(display, tabla_curenta):
    tip_de_joc = oponenti(display)
    dific = dificultate(display)

    # Pentru tipurile de joc 1 v 1 sau C vs C, tipul de estimare nu-si mai
    # are sens astfel vom omite acesta afisare in momentul alegerii setarilor jocului
    if tip_de_joc == 'cmpvspl':
        tipuri = tip_estimare(display)
    else:
        tipuri = 'nimic'

    alg = algoritm(display)

    # Daca tipul de joc este C vs C, alegerea simbolului de joc nu-si mai are rostul
    if tip_de_joc == '1v1' or tip_de_joc == 'cmpvspl':
        xand0 = x_0(display)
    else:
        xand0 = 'nimic'


    ok = Buton(display=display, top=200, left=210, w=90, h=30, text="Incepe Joc", culoareFundal=(155, 0, 55))
    ok.deseneaza()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if ok.selecteazaDupacoord(pos):
                    display.fill((0, 0, 0))  # stergere ecran
                    tabla_curenta.deseneaza_grid()
                    return xand0, alg, tip_de_joc, dific, tipuri


        pygame.display.update()


def main():


    # Vector in care memoram timpii de gandire ai calculatorului
    TIMPI_GANDIRE_PC = []


    # Timpul maxim de rulare al jocului
    TMAX = input('Introduceti TMAX:\n')
    try:
        timp_timeout = int(TMAX)
    except ValueError:
        print('Nu ati introdus un numar intreg!')
    TMAX = int(TMAX)



    # setari interf grafica
    pygame.init()
    pygame.display.set_caption("Diaconu Andrei - 10 x 10 Piese")
    # dimensiunea ferestrei in pixeli

    # Grid de 10x10
    nl = 10
    nc = 10
    w = 50
    ecran = pygame.display.set_mode(size=(nc * (w + 1) - 1, nl * (w + 1) - 1))  # N *w+ N-1= N*(w+1)-1
    Joc.initializeaza(ecran, NR_LINII=nl, NR_COLOANE=nc, dim_celula=w)


    # initializare tabla
    tabla_curenta = Joc(NR_LINII=10, NR_COLOANE=10);
    Joc.JMIN, tip_algoritm, tip_de_joc, dificultate, tipuri = incepe_joc(ecran, tabla_curenta)

    # Joc.JMIN va avea valoare 'nimic' atunci cand vom aveam C vs C
    if Joc.JMIN == 'nimic':
        Joc.JMIN = 'x'
    Joc.JMAX = '0' if Joc.JMIN == 'x' else 'x'


    print("Tabla initiala")
    print(str(tabla_curenta))


    # Stabilirea adancimii in functie de dificultate
    if dificultate == 'scazuta':
        adanc = 1
    elif dificultate == 'medie':
        adanc = 2
    elif dificultate == 'ridicata':
        adanc = 3
    stare_curenta = Stare(tabla_curenta, 'x', adancime=adanc)



    tabla_curenta.deseneaza_grid()

    # Ne va permite sa marcam prima si a doua plasare de piesa din placuta
    ok = 0

    # Vom folosi acesti indici pentru a retine linia si coloana primei plasari
    # a piesei din placuta
    linprimamutare = 0
    colprimamutare = 0

    # Stabilim timpul de inceput al jocului
    t_inainte_TMAX = int(round(time.time() * 1000))

    # Folosim o variabila tmp pentru a nu se reseta timpul de gandire al
    # utilizatorului
    tmp = 0


    nr_mutari_jucator = 0
    nr_mutari_calculator = 0


    while True:

        # Verificam daca se depaseste TMAX
         t_dupa_TMAX = int(round(time.time() * 1000))
         if t_dupa_TMAX - t_inainte_TMAX > TMAX:
             #print('S-a depasit timpul stabilit la inceputul jocului!')
             #break
             pass

        ##########################################################
        # Calculator vs jucator
         if tip_de_joc == 'cmpvspl':
             if (stare_curenta.j_curent == Joc.JMIN):
                 # if (stare_curenta.j_curent != Joc.GOL):


                 # Daca lista de mutari nu mai contine elemente insmeana ca
                 # nu se mai pot genera alte mutari si astfel se ajunge la finalul
                 # posibilitatii plasarii de simboluri
                 if len(stare_curenta.mutari()) == 0:
                     if Joc.scor_JMIN > Joc.scor_JMAX:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE

                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMIN:
                                 # Vom colora casutele castigatorului
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)

                         # Afisarile cerute
                         print('A castigat x')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: '+ str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)


                         break
                     elif Joc.scor_JMAX > Joc.scor_JMIN:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMAX:

                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat 0')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)


                         break
                     else:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('Remiza')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)


                         break

                 # Vom marca inceputul timpului de gandire al utilizatorului
                 if tmp == 0:
                     t_inainte_utilizator = int(round(time.time() * 1000))
                     tmp = 1


                 for event in pygame.event.get():
                     if event.type == pygame.QUIT:
                         # iesim din program
                         print('Scor jucator: ', Joc.scor_JMIN)
                         print('Scor calculator: ', Joc.scor_JMAX)
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ', sum(TIMPI_GANDIRE_PC)/len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ',min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL)/len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)
                         pygame.quit()
                         sys.exit()


                     elif event.type == pygame.MOUSEBUTTONDOWN:

                         pos = pygame.mouse.get_pos()  # coordonatele cursorului la momentul clickului

                         for np in range(len(Joc.celuleGrid)):

                             if Joc.celuleGrid[np].collidepoint(pos):
                                 linie = np // Joc.NR_LINII
                                 coloana = np % Joc.NR_COLOANE
                                 ###############################

                                 # Nu dam posibilitatea marcarii unui nou simbol peste altul si
                                 # de asemenea verificam sa fim la prima plasare a piesei (ok=0)
                                 if stare_curenta.tabla_joc.matr[linie][coloana] != Joc.GOL and ok == 0:
                                     continue

                                 # Prima plasare a piesei
                                 if ok == 0:
                                     linprimamutare = linie
                                     colprimamutare = coloana
                                     stare_curenta.tabla_joc.matr[linie][coloana] = stare_curenta.j_curent
                                     stare_curenta.tabla_joc.penultima_mutare = (linie, coloana)
                                     ok += 1

                                 # A doua plasare a piesei
                                 elif ok == 1:
                                     # Timpul de finalizare a plasarii placutei pentru jucator
                                     t_dupa_utilizator = int(round(time.time() * 1000))
                                     print('Utilizatorul a gandit: '+ str(t_dupa_utilizator - t_inainte_utilizator) + " milisecunde.\n")

                                     # Permitem ca jucatorul sa se razgandeasca in privinta plasarii primei piese
                                     if linie == linprimamutare and coloana == colprimamutare:
                                         stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL
                                         ok = 0

                                     # Daca vrem sa plasam al doilea simbol peste un simbol existent, nu permitem
                                     elif stare_curenta.tabla_joc.matr[linie][coloana] != Joc.GOL:
                                         continue

                                     # Daca placuta formata nu este de tipul 1x2 sau 2x1 stergem modificarile
                                     # aduse de utilizator pentru mutarea curenta si acesta are posibilitatea unei noi plasari
                                     if not ((abs(linie - linprimamutare) == 1 and coloana == colprimamutare) or (
                                             abs(coloana - colprimamutare) == 1 and linie == linprimamutare)):
                                         ok = 0
                                         stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                         stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL

                                     # Daca este o mutare valida marcam acest lucru
                                     else:
                                         stare_curenta.tabla_joc.matr[linie][coloana] = stare_curenta.j_curent
                                         stare_curenta.tabla_joc.ultima_mutare = (linie, coloana)
                                         ok += 1

                                     # Verificam acum daca intre vecinii placutei exista macr un 0 si macar un x

                                     # Consideram ca sunt pe aceasi coloana pentru a stabili natura piesei (2x1)
                                     if coloana == colprimamutare:
                                         nrx = 0
                                         nr0 = 0

                                         # Stabilim un reper de verificare a vecinilor prin selecatrea minimului dupa linie
                                         minim = min(linie, linprimamutare)
                                         dir = [[-1, 0], [-1, -1], [0, -1], [1, -1], [2, -1], [2, 0], [2, 1], [1, 1],
                                                [0, 1], [-1, 1]]

                                         # Contorizam in functie de tip simbolurile aflate in jurul placutei
                                         for item in dir:
                                             # Ne asiguram sa nu excedem dimensiunile matricei
                                             try:
                                                 if stare_curenta.tabla_joc.matr[minim + item[0]][coloana + item[1]] == Joc.JMIN:
                                                     nrx += 1
                                                 elif stare_curenta.tabla_joc.matr[minim + item[0]][
                                                     coloana + item[1]] == Joc.JMAX:
                                                     nr0 += 1
                                             except IndexError:
                                                 continue

                                         # Daca nu exista minim un 0 cat si un x nu este o mutare valida
                                         if not (nrx >= 1 and nr0 >= 1):
                                             ok = 0
                                             stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                             stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL

                                    # Procedam similar si in cazul aceleiasi linii (placuta de 1x2)
                                     elif linie == linprimamutare:
                                         nrx = 0
                                         nr0 = 0
                                         minim = min(coloana, colprimamutare)
                                         dir = [[0, -1], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [0, 2], [1, 2], [1, 1],
                                                [1, 0], [1, -1]]

                                         for item in dir:
                                             try:
                                                 if stare_curenta.tabla_joc.matr[linie + item[0]][minim + item[1]] == Joc.JMIN:
                                                     nrx += 1
                                                 elif stare_curenta.tabla_joc.matr[linie + item[0]][
                                                     minim + item[1]] == Joc.JMAX:
                                                     nr0 += 1
                                             except IndexError:
                                                 continue


                                         if not (nrx >= 1 and nr0 >= 1):
                                             ok = 0
                                             stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                             stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL

                                 # afisarea starii jocului in urma mutarii utilizatorului
                                 print("\nTabla dupa mutarea jucatorului")
                                 print(str(stare_curenta))

                                 stare_curenta.tabla_joc.deseneaza_grid()

                                 # Verificam daca ok = 2 astfel verificam corectitdinea mutarii jucatorului
                                 # si stabilim scorul
                                 if ok == 2:
                                     nr_mutari_jucator += 1
                                     if stare_curenta.j_curent == Joc.JMIN:
                                         if (afis_daca_final_pc(stare_curenta)):
                                             Joc.scor_JMIN += 1
                                         if (afis_daca_final(stare_curenta)):
                                             Joc.scor_JMIN += 1
                                         print('Scor JMIN: ', Joc.scor_JMIN)


                                 # S-a realizat o mutare. Schimb jucatorul cu cel opus
                                 if ok == 2:
                                     tmp = 0
                                     stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)
                                     ok = 0


             # --------------------------------
             # if (stare_curenta.j_curent == Joc.GOL):
             elif(stare_curenta.j_curent == Joc.JMAX):

                 if len(stare_curenta.mutari()) == 0:
                     if Joc.scor_JMIN > Joc.scor_JMAX:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMIN:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat x')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)

                         break
                     elif Joc.scor_JMAX > Joc.scor_JMIN:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMAX:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat 0')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)
                         break
                     else:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('Remiza')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar mutari jucator: ', nr_mutari_jucator)
                         print('Numar mutari calculator: ', nr_mutari_calculator)
                         break



                 # jucatorul e JMAX (calculatorul):
                 # Mutare calculator



                 # preiau timpul in milisecunde de dinainte de mutare
                 t_inainte = int(round(time.time() * 1000))



                 if tip_algoritm == 'minimax':
                     # Apelez min_max pentru starea curenta un raport cu ultima_mutare
                     stare_actualizata = min_max(stare_curenta, estimare_scor_tip=tipuri)

                     # Schimb penultima mutare ca fiind ultima pentru ai permite algoritmului sa verifice daca aduna punct si pentru aceasta
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = min_max(stare_actualizata, estimare_scor_tip=tipuri)

                     # Apelez min_max_generator pentru a genera urmatoarea mutare a calculatorului in situatiile
                     # in care functia min_max iese pe return pentru ambele piese din placuta
                     stare_actualizata = min_max_generator(stare_actualizata, estimare_scor_tip=tipuri)
                 else:  # tip_algoritm=="alphabeta"
                     # Procedam asemenator si in cazul algortimului alpha-beta
                     stare_actualizata = alpha_beta(-500, 500, stare_curenta, estimare_scor_tip=tipuri)
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = alpha_beta(-500, 500, stare_actualizata, estimare_scor_tip=tipuri)
                     stare_actualizata = alpha_beta_generator(-500, 500, stare_actualizata, estimare_scor_tip=tipuri)


                 stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc

                 # Estimarea radacinii arborelui
                 print('Scor stare curenta: ', stare_curenta.scor)

                 # Afisam nodurile generate
                 print('Numar noduri generate: ', nr_noduri_generate)

                 print("Tabla dupa mutarea calculatorului\n" + str(stare_curenta))

                 # preiau timpul in milisecunde de dupa mutare
                 t_dupa = int(round(time.time() * 1000))
                 print("Calculatorul a \"gandit\" timp de " + str(t_dupa - t_inainte) + " milisecunde.")

                 # timpii de gandire ai calculatorului
                 TIMPI_GANDIRE_PC.append(t_dupa - t_inainte)

                 stare_curenta.tabla_joc.deseneaza_grid()

                 # Stabilim scorul
                 if (afis_daca_final_pc(stare_curenta)):
                     Joc.scor_JMAX += 1
                 if (afis_daca_final(stare_curenta)):
                     Joc.scor_JMAX += 1
                 print('Scor JMAX: ', Joc.scor_JMAX)


                 nr_mutari_calculator += 1

                 # S-a realizat o mutare. Schimb jucatorul cu cel opus
                 stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)


         # Jucator vs jucator
         elif tip_de_joc == '1v1':
             #if (stare_curenta.j_curent == Joc.JMIN):
             if (stare_curenta.j_curent != Joc.GOL):

                 if len(stare_curenta.mutari()) == 0:
                     if Joc.scor_JMIN > Joc.scor_JMAX:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMIN:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat x')
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari jucatori: ', nr_mutari_jucator)
                         break


                     elif Joc.scor_JMAX > Joc.scor_JMIN:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMAX:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat 0')
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari jucatori: ', nr_mutari_jucator)
                         break


                     else:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('Remiza')
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari jucatori: ', nr_mutari_jucator)
                         break


                 if tmp == 0:
                     t_inainte_utilizator = int(round(time.time() * 1000))
                     tmp = 1


                 for event in pygame.event.get():
                     if event.type == pygame.QUIT:
                         t_final_joc = int(round(time.time() * 1000))
                         print('Scor jucator1: ', Joc.scor_JMIN)
                         print('Scor jucator2: ', Joc.scor_JMAX)
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari jucatori: ', nr_mutari_jucator)
                         # iesim din program
                         pygame.quit()
                         sys.exit()


                     elif event.type == pygame.MOUSEBUTTONDOWN:

                         pos = pygame.mouse.get_pos()  # coordonatele cursorului la momentul clickului

                         for np in range(len(Joc.celuleGrid)):

                             if Joc.celuleGrid[np].collidepoint(pos):
                                 linie = np // Joc.NR_LINII
                                 coloana = np % Joc.NR_COLOANE
                                 ###############################


                                 if stare_curenta.tabla_joc.matr[linie][coloana] != Joc.GOL and ok == 0:
                                     continue


                                 if ok == 0:
                                     linprimamutare = linie
                                     colprimamutare = coloana
                                     stare_curenta.tabla_joc.matr[linie][coloana] = stare_curenta.j_curent
                                     stare_curenta.tabla_joc.ultima_mutare = (linie, coloana)
                                     stare_intermediara1 = copy.deepcopy(stare_curenta)
                                     ok += 1


                                 elif ok == 1:
                                     t_dupa_utilizator = int(round(time.time() * 1000))
                                     print('Utilizatorul a gandit: ' + str(t_dupa_utilizator - t_inainte_utilizator) + " milisecunde.\n")

                                     if linie == linprimamutare and coloana == colprimamutare:
                                         stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL
                                         ok = 0
                                     elif stare_curenta.tabla_joc.matr[linie][coloana] != Joc.GOL:
                                         continue

                                     if not ((abs(linie - linprimamutare) == 1 and coloana == colprimamutare) or (
                                             abs(coloana - colprimamutare) == 1 and linie == linprimamutare)):
                                         ok = 0
                                         stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                         stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL
                                     else:
                                         stare_curenta.tabla_joc.matr[linie][coloana] = stare_curenta.j_curent
                                         stare_curenta.tabla_joc.ultima_mutare = (linie, coloana)
                                         stare_intermediara2 = copy.deepcopy(stare_curenta)
                                         ok += 1

                                     if coloana == colprimamutare:
                                         nrx = 0
                                         nr0 = 0
                                         minim = min(linie, linprimamutare)
                                         dir = [[-1, 0], [-1, -1], [0, -1], [1, -1], [2, -1], [2, 0], [2, 1], [1, 1],
                                                [0, 1], [-1, 1]]
                                         for item in dir:
                                             try:
                                                 if stare_curenta.tabla_joc.matr[minim + item[0]][coloana + item[1]] == Joc.JMIN:
                                                     nrx += 1
                                                 elif stare_curenta.tabla_joc.matr[minim + item[0]][
                                                     coloana + item[1]] == Joc.JMAX:
                                                     nr0 += 1
                                             except IndexError:
                                                 continue
                                         if not (nrx >= 1 and nr0 >= 1):
                                             ok = 0
                                             stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                             stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL

                                     elif linie == linprimamutare:
                                         nrx = 0
                                         nr0 = 0
                                         minim = min(coloana, colprimamutare)
                                         dir = [[0, -1], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [0, 2], [1, 2], [1, 1],
                                                [1, 0], [1, -1]]
                                         for item in dir:
                                             try:
                                                 if stare_curenta.tabla_joc.matr[linie + item[0]][minim + item[1]] == Joc.JMIN:
                                                     nrx += 1
                                                 elif stare_curenta.tabla_joc.matr[linie + item[0]][
                                                     minim + item[1]] == Joc.JMAX:
                                                     nr0 += 1
                                             except IndexError:
                                                 continue
                                         if not (nrx >= 1 and nr0 >= 1):
                                             ok = 0
                                             stare_curenta.tabla_joc.matr[linprimamutare][colprimamutare] = Joc.GOL
                                             stare_curenta.tabla_joc.matr[linie][coloana] = Joc.GOL

                                 # afisarea starii jocului in urma mutarii utilizatorului
                                 print("\nTabla dupa mutarea jucatorului")
                                 print(str(stare_curenta))

                                 stare_curenta.tabla_joc.deseneaza_grid()

                                 if ok == 2:
                                     if stare_curenta.j_curent == Joc.JMIN:
                                         if (afis_daca_final(stare_intermediara1)):
                                             Joc.scor_JMIN += 1
                                         if (afis_daca_final(stare_intermediara2)):
                                             Joc.scor_JMIN += 1
                                         print('Scor JMIN: ', Joc.scor_JMIN)
                                     elif stare_curenta.j_curent == Joc.JMAX:
                                         if (afis_daca_final(stare_intermediara1)):
                                             Joc.scor_JMAX += 1
                                         if (afis_daca_final(stare_intermediara2)):
                                             Joc.scor_JMAX += 1
                                         print('Scor JMAX: ', Joc.scor_JMAX)


                                 # S-a realizat o mutare. Schimb jucatorul cu cel opus
                                 if ok == 2:
                                     nr_mutari_jucator += 1
                                     tmp = 0
                                     stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)
                                     ok = 0

         # Calculator vs calculator
         elif tip_de_joc=='pc1vspc2':
             if (stare_curenta.j_curent == Joc.JMIN):
                 if len(stare_curenta.mutari()) == 0:
                     if Joc.scor_JMIN > Joc.scor_JMAX:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMIN:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat x')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)


                         break
                     elif Joc.scor_JMAX > Joc.scor_JMIN:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMAX:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat 0')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)


                         break
                     else:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('Remiza')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)


                         break

                 # jucatorul e JMAX (calculatorul):
                 # Mutare calculator

                 # preiau timpul in milisecunde de dinainte de mutare
                 t_inainte = int(round(time.time() * 1000))
                 if tip_algoritm == 'minimax':
                     stare_actualizata = min_max(stare_curenta, estimare_scor_tip='tip1')
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = min_max(stare_actualizata, estimare_scor_tip='tip1')
                     stare_actualizata = min_max_generator(stare_actualizata, estimare_scor_tip='tip1')

                 else:  # tip_algoritm=="alphabeta"
                     stare_actualizata = alpha_beta(-500, 500, stare_curenta, estimare_scor_tip='tip1')
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = alpha_beta(-500, 500, stare_actualizata, estimare_scor_tip='tip1')
                     stare_actualizata = alpha_beta_generator(-500, 500, stare_actualizata, estimare_scor_tip='tip1')


                 stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc
                 print('Numar noduri generate: ', nr_noduri_generate)

                 # Estimarea radacinii arborelui
                 print('Scor stare curenta: ', stare_curenta.scor)

                 print("Tabla dupa mutarea calculatorului\n" + str(stare_curenta))

                 # preiau timpul in milisecunde de dupa mutare
                 t_dupa = int(round(time.time() * 1000))
                 print("Calculatorul a \"gandit\" timp de " + str(t_dupa - t_inainte) + " milisecunde.")


                 TIMPI_GANDIRE_PC.append(t_dupa - t_inainte)


                 stare_curenta.tabla_joc.deseneaza_grid()
                 if stare_curenta.j_curent == Joc.JMIN:
                     if (afis_daca_final_pc(stare_curenta)):
                         Joc.scor_JMIN += 1
                     if (afis_daca_final(stare_curenta)):
                         Joc.scor_JMIN += 1
                     print('Scor JMIN: ', Joc.scor_JMIN)

                 # S-a realizat o mutare. Schimb jucatorul cu cel opus
                 nr_mutari_calculator += 1
                 stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)

             else:
                 if len(stare_curenta.mutari()) == 0:
                     if Joc.scor_JMIN > Joc.scor_JMAX:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMIN:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat x')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)
                         break


                     elif Joc.scor_JMAX > Joc.scor_JMIN:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             if stare_curenta.tabla_joc.matr[i][j] == Joc.JMAX:
                                 stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)


                         print('A castigat 0')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)
                         break


                     else:
                         for item in range(len(Joc.celuleGrid)):
                             i = item // Joc.NR_LINII
                             j = item % Joc.NR_COLOANE
                             stare_curenta.tabla_joc.deseneaza_castigator(lin=i, col=j)
                         print('Remiza')
                         print('Timpul minim de gandire al calculatorului: ', min(TIMPI_GANDIRE_PC))
                         print('Timpul maxim de gandire al calculatorului: ', max(TIMPI_GANDIRE_PC))
                         print('Timpul mediu de gandire al calculatorului: ',sum(TIMPI_GANDIRE_PC) / len(TIMPI_GANDIRE_PC))
                         print('Mediana corespunzatoare timpilor de gandire: ', median(TIMPI_GANDIRE_PC))
                         print('Numar minim noduri generate: ', min(VECT_NR_TOTAL))
                         print('Numar maxim noduri generate: ', max(VECT_NR_TOTAL))
                         print('Numar mediu noduri generate: ', sum(VECT_NR_TOTAL) / len(VECT_NR_TOTAL))
                         print('Mediana noduri generate: ', median(VECT_NR_TOTAL))
                         t_final_joc = int(round(time.time() * 1000))
                         print('Timp fial joc: ' + str(t_final_joc - t_inainte_TMAX))
                         print('Numar total mutari: ', nr_mutari_calculator)
                         break

                 # jucatorul e JMAX (calculatorul):
                 # Mutare calculator

                 # preiau timpul in milisecunde de dinainte de mutare
                 t_inainte = int(round(time.time() * 1000))
                 if tip_algoritm == 'minimax':
                     stare_actualizata = min_max(stare_curenta, estimare_scor_tip='tip2')
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = min_max(stare_actualizata, estimare_scor_tip='tip2')
                     stare_actualizata = min_max_generator(stare_actualizata, estimare_scor_tip='tip2')


                 else:  # tip_algoritm=="alphabeta"
                     stare_actualizata = alpha_beta(-500, 500, stare_curenta, estimare_scor_tip='tip2')
                     stare_actualizata.tabla_joc.ultima_mutare = stare_actualizata.tabla_joc.penultima_mutare
                     stare_actualizata = alpha_beta(-500, 500, stare_actualizata, estimare_scor_tip='tip2')
                     stare_actualizata = alpha_beta_generator(-500, 500, stare_actualizata, estimare_scor_tip='tip2')


                 stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc

                 # Estimarea radacinii arborelui
                 print('Scor stare curenta: ', stare_curenta.scor)

                 print('Numar noduri generate: ', nr_noduri_generate)

                 print("Tabla dupa mutarea calculatorului\n" + str(stare_curenta))

                 # preiau timpul in milisecunde de dupa mutare
                 t_dupa = int(round(time.time() * 1000))
                 print("Calculatorul a \"gandit\" timp de " + str(t_dupa - t_inainte) + " milisecunde.")


                 TIMPI_GANDIRE_PC.append(t_dupa - t_inainte)


                 stare_curenta.tabla_joc.deseneaza_grid()
                 if stare_curenta.j_curent == Joc.JMAX:
                     if (afis_daca_final_pc(stare_curenta)):
                         Joc.scor_JMAX += 1
                     if (afis_daca_final(stare_curenta)):
                         Joc.scor_JMAX += 1
                     print('Scor JMAX: ', Joc.scor_JMAX)

                 # S-a realizat o mutare. Schimb jucatorul cu cel opus
                 nr_mutari_calculator += 1
                 stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)

if __name__ == "__main__":
    main()
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Autor: Marian Kravec (2mAIN)

# Stabilné parametre:
n = 1000 # pocet simulacii
delta = 600 # cas prijmania novych zakaznikov
mean = 1 # stredna hodnota gamma rozdelenie
H = np.log(2)/10 # parameter pre pravdepodobnost vstupu do fronty
cl = 0.95 # confidence level intervalu spolahlivosti

# Nastavované parametre:
lambdas = [0.5, 1, 2] # mozne hodnoty lambda
var_coefs = [0.5, 1, 2] # mozne hodnoty koeficientu variacie

fig, axs = plt.subplots(3,3)
fig.suptitle('Histogramy počtu obslúžených zákazníkov')

for l, lambd in enumerate(lambdas):
    for v, var_coef in enumerate(var_coefs):

        # Pomocné výpočty :
        # Výpočet parametrov pre Gamma rozdelenie

        # var_coef = std/mean
        std = mean * var_coef

        # mean = shape * scale
        # std = sqrt(shape) * scale

        # mean/std = (shape * scale)/(sqrt(shape) * scale)
        #          = shape/sqrt(shape)
        #          = sqrt(shape)
        shape = (mean/std)**2

        # mean/(std^2) = (shape * scale)/(shape * scale^2)
        #              = scale/(scale^2)
        #              = scale^(-1)
        scale = (mean/(std**2))**(-1)

        # Inicializacia zoznamov výsledkov jednotlivých runov
        zoz_obsluzeny = []
        zoz_max_fronta = []
        zoz_priem_cas = []

        # Simulacia
        for i in range(n):
            # vynulovanie hodnot na zaciatku dna
            fronta = 0 # pocet ludi vo fronte
            obsluzeny = 0 # pocet obsluzenych
            sucasny_cas = 0 # sucasny cas
            cas_posl_obs = 0 # cas kedy bol obsluzeny posledny zakaznik
            max_fronta = 0 # dlzka maximalnej fronty
            celkova_fronta = 0 # sucet vsetkych dlzok fronty (pre vypocet priemernej)


            zakaznik_prich = np.random.gamma(shape, scale, 1)[0] # prichod prveho zakaznika
            while sucasny_cas + zakaznik_prich <= delta: # ak zakaznik prisiel v akceptacnom case
                sucasny_cas += zakaznik_prich # posunitie sucasneho casu o prichod zakaznika
                while sucasny_cas > cas_posl_obs+lambd: # ak za ten cas stihli inych zakaznikov obsluzit
                    if fronta == 0: # ak nikto nie je vo fronte
                        cas_posl_obs = sucasny_cas # tak zakaznik moze byt obsluzeny najskor o lambda
                                                   # preto posunieme umelo cas posledneho obsluzeneho na sucasny cas
                    else: # ak je niekto vo fronte
                        cas_posl_obs += lambd # posunieme cas posledneho obsluzeneho
                        fronta -= 1 # odoberieme z fronty
                        obsluzeny += 1 # pridame k obsluzenym

                prav_do_fronty = np.exp(-H*fronta) # vypocet pravdepodobnosti vstupu do fronty
                rozhod_zakaz = np.random.uniform(0, 1, 1)[0] # rozhodnutie zakaznika
                if rozhod_zakaz <= prav_do_fronty: # ak je rozhodnutie menej ako pravdepodobnost vstupy do fronty
                    fronta += 1
                    celkova_fronta += fronta

                max_fronta = max(max_fronta, fronta) # update maximalnej fronty na zaklade sucasnej
                zakaznik_prich = np.random.gamma(shape, scale, 1)[0] # prichod dalsieho zakaznika

            obsluzeny += fronta # pridanie obsluzenie ludi vo fronte

            priem_fronta = celkova_fronta/obsluzeny # priemerna dlzka fronty ako celkova dlzka fronty ked sa zakaznici zaradili
                                                    # vydelena poctom obsluzenych zakazniko
            priem_cas = priem_fronta*lambd # priemerny cas straveny vo fronte ako priemrna dlzky fronty vynasobena dlzkov cakania a vyriešenie jedneho

            #print(obsluzeny, max_fronta, np.round(priem_cas, 2))

            #sucty vsetkych dni
            zoz_obsluzeny.append(obsluzeny)
            zoz_max_fronta.append(max_fronta)
            zoz_priem_cas.append(priem_cas)

        # vypocty priemernych hodnot:
        oblsuzeny_priem = np.sum(zoz_obsluzeny)/n
        max_fronta_priem = np.sum(zoz_max_fronta)/n
        priem_cas_priem = np.sum(zoz_priem_cas)/n

        # standardna odchylka
        oblsuzeny_sem = stats.sem(zoz_obsluzeny)
        max_fronta_sem = stats.sem(zoz_max_fronta)
        priem_cas_sem = stats.sem(zoz_priem_cas)

        # vypocty cl percentneho intervalu spolahlivosti
        oblsuzeny_conf_int = np.round(stats.t.interval(confidence=cl, df=n-1, loc=oblsuzeny_priem, scale=oblsuzeny_sem), 3)
        max_fronta_conf_int = np.round(stats.t.interval(confidence=cl, df=n-1, loc=max_fronta_priem, scale=max_fronta_sem), 3)
        priem_cas_conf_int = np.round(stats.t.interval(confidence=cl, df=n-1, loc=priem_cas_priem, scale=priem_cas_sem), 3)

        # vypisy vysledkou
        print("Lambda: {}, Koeficient variacie: {}".format(lambd, var_coef))
        print("Informácia                         Priemer  Interval spolahlivosti  Šírka intervalu")
        print("Priemerny pocet obsluzench:       {:8.3f}    {:20}  {}".format(oblsuzeny_priem, str(oblsuzeny_conf_int), np.round(oblsuzeny_conf_int[1]-oblsuzeny_conf_int[0], 3)))
        print("Priemrna dlzka maximalnej fronty: {:8.3f}    {:20}  {}".format(max_fronta_priem, str(max_fronta_conf_int), np.round(max_fronta_conf_int[1]-max_fronta_conf_int[0], 3)))
        print("Priemerny cas cakanie:            {:8.3f}    {:20}  {}".format(priem_cas_priem, str(priem_cas_conf_int), np.round(priem_cas_conf_int[1]-priem_cas_conf_int[0], 3)))


        print()

        axs[l, v].title.set_text('Lamda: {} Koef. var.:{}'.format(lambd, var_coef))
        #axs[l, v].set_xlim([250, 650])
        axs[l, v].hist(zoz_obsluzeny)

fig.tight_layout()
plt.show()

"""
Komentar k vysledkom:

Vidime, ze zo zvysujucou sa labdou sa nam znisuje pocet obsluzenych a zvysuje dlzka fronty 
a dlzka cakanie co je ocakavatelne kedze cakanie jedneho zakaznika trva dlhsie.
Zaroven si ale mozeme vsimnut, ze intervali spolahlivosti sa spravaju zvlastne,
interval poctu zakaznikov sa zmesuje, zatial co intervaly dlzky fronty a casu sa zvacsuju 
co je asi sposobene tym, ze sa samotna hodnota poctu zakaznikov znizuje a dlzky fronty a case
zvysuje... tu si nie som isty.

V pripade koeficientu variacie si mozeme vsimnut, ze z jeho zvysujucou hodnotou klesa priemerny
pocet obsluzenych a rastie dlzka fronty a cas cakania, toto je s najvacsou pravdepodobnostou
sposobene tym, ze so zvysujucou sa variaciou sa castejsie stava ze pride viac zakaznikov v kratkom
case co predlzi frontu a cas a zaroven zvysi sancu ze zakaznik do frontu nevstupy, pricom tento
prichod vela zakaznikov sa strieda o obdobim kedy nikto nepride. Ak ide o intervaly spolahlivosti tak
tie sa zo zvacsujucou variaciou ocakavatelne zvacsuju
"""
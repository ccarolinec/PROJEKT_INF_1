**________TRANFORMACJE WSPÓŁRZĘDNYCH________**

Projekt: Informatyka Geodezyjna II

Podany program służy do tranformacji współrzędnych między układami geocentrycznym XYZ, elipsoidalnym (phi, lambda h), topocentrycznym NEU oraz układami 1992 i 2000 na elipsoidach GRS80 i WGS84.

Wymagania:
- program python w wersji 3.11
- biblioteka math
- biblioteka numpy 
- biblioteka sys
- system operacyjny Windows 11

__Instrukcja użycia programu:__

1. Program jest obsługiwany w Terminalu i działa za pomocą podawania konkretnych flag. Każde użycie zaczyna się od wpisania python i wpisaniu pliku prgramu - geo_v1.py:

        python geo_v1.py
  
2. Następnie należy uwzględnić liczbę linijek nagłówka flagą --header_line w pliku wejściowym ze współrzędnymi tak, aby progam zaczął sczytywanie od danej linijki. Podajemy liczbę wersów, które zajmuje nagłówek:

        python geo_v1.py --header_lines 4

3. Potem uzytkownik podaje flagę funkcji obliczeniowej, którą zamierzamy wykonać:

- transformacja ze współrzędnych X,Y,Z do współrzędnych phi, lambda, wysokość:

      python geo_v1.py --header_lines 4 --xyz2plh

- transformacja ze współrzędnych phi, lambda, wysokość do wspołrzędnych X,Y,Z:

        python geo_v1.py --header_lines 4 --plh2xyz

- transformacja ze współrzędnych X,Y,Z do układu N,E,U (w tym przypadku należy także podać współrzędne geocentryczne anteny)

      python geo_v1.py --header_lines 4 --xyz2neu 3664940.500 1409153.590 5009571.170

- transformacja ze współrzędnych phi, lamda do układu 2000

        python geo_v1.py --header_lines 4 --pl22000

- transformacja ze współrzędnych phi, lambda do układu 1992

      python geo_v1.py --header_lines 4 --pl21992

4. Nastepnie użytkownik jest proszony o podanie powierzchni odniesienia - do wyboru są elipsoida Krasowskiego, WRS84 i GRS80. Program przyjmuje jedynie pisownię małymi literami.

           python geo_v1.py --header_lines 4 --xyz2plh wrs84

           python geo_v1.py --header_lines 4 --xyz2plh grs80
   
           python geo_v1.py --header_lines 4 --xyz2plh Krasowski 

6. Na końcu uwzględniamy plik wejściowy o formacie .txt ze współrzędnymi. Plik ten powinien znajdować się na komputerze w folderze z plikiem geo_v1.py.
   
         python geo_v1.py --header_lines 4 --xyz2plh grs80 wsp_inp.txt
   
Format pliku: 
   W przypadku transformacji xyz2plh i xyz2neu, w w pierwszej kolumnie ma znajdować się współrzędna X, w drugiej Y, w trzeciej Z wyrażone w metrach. 
   Przy transformacji plh2xyz w trzech kolejnych kolumnach znajdują się phi, lambda i wysokość. Współrzędne te powinny być stopniach dziesiętnych.
   Podczas transformacji pl22000 i pl21992 program pobiera jedynie współrzędne phi i lambda w kolejnych kolumach wyrażone w stopniach dziesiętnych, ale plik wejściowy może zawierać informacje o wysokości - program je zignoruje.

Ułamkowe części współrzędnych wpisujemy po kropce.

6. Po kliknięciu ENTER i odświeżeniu danego folderu powinien pojawić się plik result.txt, gdzie w nazwie pliku jest uwzględniony rodzaj transformacji. Jeśli wynikiem są współrzędne wyrażone w metrach, program zaokrągli je do trzech miejsc po przecinku, czyli do milimetów.

**!!!**
Istnieje także możliwość, aby program wykonał dwie transformacje na raz: xyz2plh i xyz2neu oraz pl21990 i pl2200. 

    python geo_v1.py --header_lines 4 --xyz2plh --xyz2neu 3664940.500 1409153.590 5009571.170

    python geo_v1.py --header_lines 4 --pl21992 --pl22000 

Istotna jest kolejność wpisywania transformacji, kiedy wpiszemy je na odwrót, program wykona tylko drugą w kolejności transformację. Wyniki podwójnej transformacji zapiszą się w jednym pliku .txt.

# MashGraph

В main.cpp реализована программа трассировки лучей. Выполнены примитивы: куб, бескочесная плоскость с текстурой и сфера. Материалы объектов могут отражать и преломлять. Для уменьшения времени работы программы используется многопоточность.

Порядок компиляции:

mkdir build
cd build
cmake −DCMAKE_BUILD_TYPE=Release ..
make −j 4

./rt −out <output_path> −scene 1 −threads <threads>

где output_path - путь к выходному изображению (относительный).
 threads - количество потоков.

На выходе получаем bmp изображение (есть возможность выводить в ppm формате) размером 1280х720 пикселей. 

Используемая литература:
https://habr.com/en/post/333932/
https://habr.com/en/post/353054/
https://habr.com/ru/post/436790/
https://en.wikipedia.org/wiki/Phong_reflection_model


------------------------------------------------------------------------
The implementation of the foundations of three-dimensional mathematics and the basic concepts of computer graphics. Used Books:

1. https://habr.com/en/post/333932/
2. https://habr.com/en/post/353054/
3. https://habr.com/ru/post/436790/
4. https://en.wikipedia.org/wiki/Phong_reflection_model


#include <iostream>
#include <cmath>
#include <vector>
#include "SFML/Graphics.hpp"
const int amountPoints = 20;
const int numberLayers = 2;
const double G = 6.67408 * pow(10, -11);
const double earthMass = 6 * pow(10, 24);
const double earthMassPoint = earthMass / (amountPoints * numberLayers);
const double moonMass = 7.36 * pow(10, 22);
const double K = 0.0000000001;
const double h = 10000;

struct Coordinate{
    Coordinate(double x, double y){
        this->x = x;
        this->y = y;
    }
    Coordinate():x(0.0), y(0.0){}
    double x;
    double y;
    friend Coordinate operator * (const double & digit, const Coordinate & right);
    friend Coordinate operator * (const Coordinate & left, const Coordinate & right);
    friend Coordinate operator + (const Coordinate & left, const Coordinate & right);
    friend Coordinate operator + (const double & digit, const Coordinate &right);
    friend std::ostream & operator << (std::ostream & os, Coordinate & position);
};

std::ostream &operator<<(std::ostream &os, Coordinate &position) {
    os<<"(x, y) "<<position.x<<" "<<position.y<<"\n";
    return os;
}

Coordinate operator*(const double & digit, const Coordinate &right) {
    Coordinate answer;
    answer.x = right.x * digit;
    answer.y = right.y * digit;
    return answer;
}

Coordinate operator+(const double & digit, const Coordinate &right) {
    Coordinate answer;
    answer.x = right.x + digit;
    answer.y = right.y + digit;
    return answer;
}

Coordinate operator+(const Coordinate& left, const Coordinate& right) {
    Coordinate answer;
    answer.x = left.x + right.x;
    answer.y = left.y + right.y;
    return answer;
}

Coordinate operator*(const Coordinate &left, const Coordinate &right) {
    Coordinate answer;
    answer.x = left.x * right.x;
    answer.y = left.y * right.y;
    return answer;
}

Coordinate operator - (const Coordinate& a, const Coordinate& b){
    return {a.x - b.x, a.y - b.y};
}

double ComputeDistance(const Coordinate & pos1,const Coordinate & pos2){
    return sqrt(pow(pos2.x - pos1.x, 2) + pow(pos2.y - pos1.y, 2));
}

Coordinate ComputeElasticity(std::vector<Coordinate> & earthPoints, int numberCurrentPoint, std::vector<Coordinate> & originalPositionPoints,
                                          int numberNeighbor){
    Coordinate answer, vectorDistance;
    double distance, originalDistance;

    vectorDistance = earthPoints[numberNeighbor] - earthPoints[numberCurrentPoint];
    distance = ComputeDistance(earthPoints[numberNeighbor], earthPoints[numberCurrentPoint]);
    originalDistance = ComputeDistance(originalPositionPoints[numberNeighbor],
                                       originalPositionPoints[numberCurrentPoint]);
    answer = (distance - originalDistance) * (K / distance) * vectorDistance;

    return answer;
}



struct States{
    std::vector<Coordinate> positionEarthPoints;
    std::vector<Coordinate> speedEarthPoints;
    std::vector<Coordinate> positionEarthPointsOriginal;
    Coordinate positionMoon;
    Coordinate speedMoon;
    States(std::vector<Coordinate> & newPositionFirst, std::vector<Coordinate> & newSpeedFirst, Coordinate & newPositionSecond, Coordinate & newSpeedSecond){
        // first point will be center of earth
        positionEarthPoints = newPositionFirst;
        positionEarthPointsOriginal = newPositionFirst;
        speedEarthPoints = newSpeedFirst;
        positionMoon = newPositionSecond;
        speedMoon = newSpeedSecond;
    }
    States(const States & other){
        this->positionEarthPointsOriginal = other.positionEarthPointsOriginal;
        this->positionEarthPoints = other.positionEarthPoints;
        this->positionMoon = other.positionMoon;
        this->speedEarthPoints = other.speedEarthPoints;
        this->speedMoon = other.speedMoon;
    }
    States() = default;

    Coordinate Elasticity(std::vector<Coordinate> & earthPoints, int numberCurrentPoint, std::vector<Coordinate> & originalPositionPoints,
                          Coordinate & earthCentre){
        Coordinate answer, vectorDistance;
        double distance, originalDistance;
        int numberNeighborPoint;
        // right neighbor
        if (numberCurrentPoint % amountPoints < amountPoints - 1) {
            numberNeighborPoint = numberCurrentPoint + 1;
        } else {
            numberNeighborPoint = numberCurrentPoint - (amountPoints - 1);
        }
        answer = answer +
                 ComputeElasticity(earthPoints, numberCurrentPoint, originalPositionPoints, numberNeighborPoint);
        // left neighbor
        if (numberCurrentPoint % amountPoints > 0) {
            numberNeighborPoint = numberCurrentPoint - 1;
        } else {
            numberNeighborPoint = numberCurrentPoint + amountPoints - 1;
        }
        answer = answer +
                 ComputeElasticity(earthPoints, numberCurrentPoint, originalPositionPoints, numberNeighborPoint);
        // on next layer neighbor
        if (numberCurrentPoint + amountPoints >= earthPoints.size()) {
            vectorDistance = earthCentre - earthPoints[numberCurrentPoint];
            distance = ComputeDistance(earthCentre, earthPoints[numberCurrentPoint]);
            originalDistance = ComputeDistance(this->positionEarthPointsOriginal[0],
                                               originalPositionPoints[numberCurrentPoint]);
            answer = answer + (distance - originalDistance) * (K / distance) * vectorDistance;
        } else {
            numberNeighborPoint = amountPoints + numberCurrentPoint;
            answer = answer + ComputeElasticity(earthPoints, numberCurrentPoint, originalPositionPoints,
                                                numberNeighborPoint);
        }
        // on previous layer neighbor
        if (numberCurrentPoint - amountPoints >= 0) {
            numberNeighborPoint = numberCurrentPoint - amountPoints;
            answer = answer + ComputeElasticity(earthPoints, numberCurrentPoint, originalPositionPoints,
                                                numberNeighborPoint);
        }

        return answer;
    }

    void Derivative(){
        States copy = *this;
        for(int number = 0; number < positionEarthPoints.size(); number ++){
            Coordinate vectorDistance = positionMoon - positionEarthPoints[number];
            double distance = ComputeDistance(positionEarthPoints[number], positionMoon);
            positionEarthPoints[number] = speedEarthPoints[number];
            Coordinate elastic;
            // compute for all points , except center
            if (number != 0) {
                int numberCpy = number - 1;
                std::vector<Coordinate> cpyPosition = copy.positionEarthPoints;
                std::vector<Coordinate> cpyPointsOriginal = positionEarthPointsOriginal;
                cpyPosition.erase(cpyPosition.begin());
                cpyPointsOriginal.erase(cpyPointsOriginal.begin());
                elastic = Elasticity(cpyPosition, numberCpy, cpyPointsOriginal, copy.positionEarthPoints[0]);
            }
            speedEarthPoints[number] = Coordinate(( G * moonMass * vectorDistance.x / pow(distance, 3)),
                                                  ( G * moonMass * vectorDistance.y / pow(distance, 3))) + elastic;
        }
        // compute for moon
        Coordinate vectorDistance = copy.positionEarthPoints[0] - positionMoon;
        double distance = ComputeDistance(copy.positionEarthPoints[0], positionMoon);
        positionMoon = speedMoon;
        speedMoon = Coordinate(( G * earthMass * vectorDistance.x / pow(distance, 3)),
                               ( G * earthMass * vectorDistance.y / pow(distance, 3)));
    }

    friend States operator * (double digit,const States& right);
    friend States operator + (const States&  left, const States& right);
    friend std::ostream & operator << (std::ostream & os, States & states);
};

std::vector<Coordinate> operator * (double digit, const std::vector<Coordinate>& right){
    std::vector<Coordinate> answer;
    answer.reserve(right.size());
    for(auto i : right)
        answer.push_back(digit * i);
    return answer;
}

std::vector<Coordinate> operator + (const std::vector<Coordinate>& left, const std::vector<Coordinate>& right){
    std::vector<Coordinate> answer;
    answer.reserve(right.size());
    for(int i=0; i <right.size(); i++)
        answer.push_back(left[i] + right[i]);
    return answer;
}

std::ostream &operator<<(std::ostream &os, States &states) {
    for(int i = 0; i <states.positionEarthPoints.size(); i++) {
        os<<i<< "   Position Earth - " << states.positionEarthPoints[i] << "\n"
          <<  i<< "   Speed Earth - " << states.speedEarthPoints[i] << "\n";
    }
    os << "Position Moon - " << states.positionMoon << "\n"
       << "Speed Moon - " << states.speedMoon << "\n";
    return os;
}

States operator*(double digit, const States& right) {
    States answer;
    answer.positionEarthPoints = digit * right.positionEarthPoints;
    answer.speedEarthPoints = digit * right.speedEarthPoints;
    answer.positionMoon = digit * right.positionMoon;
    answer.speedMoon = digit * right.speedMoon;
    answer.positionEarthPointsOriginal = right.positionEarthPointsOriginal;
    return answer;
}

States operator+(const States& left, const States& right) {
    States answer;
    answer.positionEarthPoints = left.positionEarthPoints + right.positionEarthPoints;
    answer.speedEarthPoints = left.speedEarthPoints + right.speedEarthPoints;
    answer.positionMoon = left.positionMoon + right.positionMoon;
    answer.speedMoon = left.speedMoon + right.speedMoon;
    answer.positionEarthPointsOriginal = right.positionEarthPointsOriginal;
    return answer;
}

void Runge_Kutta(States & inputState, double h){
    States k1;
    k1 = inputState;
    k1.Derivative();
    States k2 = inputState + h/2 * k1;
    k2.Derivative();
    States k3 = inputState + h/2 * k2;
    k3.Derivative();
    States k4 = inputState + h * k3;
    k4.Derivative();
    inputState = inputState + h/6  * (k1 + (2 * k2) + (2 * k3) +k4);
}

void Barycenter(States & inputParams) {
    double xSumMass=0, ySumMass=0, xSumSpeed=0, ySumSpeed=0;

    xSumMass = moonMass * inputParams.positionMoon.x;
    ySumMass = moonMass * inputParams.positionMoon.y;
    xSumSpeed = moonMass * inputParams.speedMoon.x;
    ySumSpeed = moonMass * inputParams.speedMoon.y;

    for(int i=0; i < inputParams.positionEarthPoints.size(); i++) {
        xSumMass = xSumMass + earthMassPoint * inputParams.positionEarthPoints[i].x;
        ySumMass = ySumMass + earthMassPoint * inputParams.positionEarthPoints[i].y;
    }
    xSumMass = xSumMass/(earthMass + moonMass);
    ySumMass = ySumMass/(earthMass + moonMass);
    xSumSpeed = xSumSpeed/(earthMass + moonMass);
    ySumSpeed = ySumSpeed/(earthMass + moonMass);

    for(int i=0; i < inputParams.positionEarthPoints.size(); i++) {
        inputParams.positionEarthPoints[i].x -= xSumMass;
        inputParams.positionEarthPoints[i].y -= ySumMass;
        inputParams.speedEarthPoints[i].x -= xSumSpeed;
        inputParams.speedEarthPoints[i].y -= ySumSpeed;
    }

    inputParams.positionMoon.x -= xSumMass;
    inputParams.positionMoon.y -= ySumMass;
    inputParams.speedMoon.x -= xSumSpeed;
    inputParams.speedMoon.y -= ySumSpeed;
}

void CreatePointEarth(std::vector<Coordinate> & vectorPoints){
    int center = 5 * pow(10, 8);
    int radius = 6371 * pow(10, 3);
    for(int j=1; j<=numberLayers; j++){
        for (int i=0; i < amountPoints; i++){
            Coordinate point ((double)(center + radius/j * cos(i * M_PI * 2 / amountPoints)), ((double)round(center + radius / j * sin(i * M_PI * 2 / amountPoints))));
            vectorPoints.push_back(point);
        }
    }
}

int main() {
    sf::RenderWindow window;
    window.create(sf::VideoMode(1000, 1000), "MoonAndEarth");
    sf::View view;
    view.setSize(2000000000, 2000000000);
    int radiusMoon = 50000000, radiusEarth = 63710000;
    sf::CircleShape Moon(radiusMoon);
    Moon.setFillColor(sf::Color::Red);
    window.setFramerateLimit(100);

    Coordinate centreEarth(500000000, 500000000);
    Coordinate positionMoon(862000000, 500000000);
    Coordinate speedMoon(0, 970);

    std::vector<Coordinate> positionPointsEarth = {centreEarth};
    CreatePointEarth(positionPointsEarth);

    std::vector<Coordinate> speedPointsEarth (positionPointsEarth.size(), Coordinate(0, 0));
    std::vector<sf::CircleShape> earthCircles;

    for(int i=0; i < positionPointsEarth.size(); i++){
        if (i == 0) {
            earthCircles.emplace_back(radiusEarth/10);
            earthCircles[i].setFillColor(sf::Color::Green);
        }
        else{
            earthCircles.emplace_back(radiusEarth/10);
            earthCircles[i].setFillColor(sf::Color::Blue);
        }
    }

    States stateStartParams(positionPointsEarth, speedPointsEarth, positionMoon, speedMoon);
    Barycenter(stateStartParams);
    view.setCenter(stateStartParams.positionEarthPoints[0].x, stateStartParams.positionEarthPoints[0].y);
    int N = 1000000;
    window.setView(view);
    while (N > 0) {
        Runge_Kutta(stateStartParams, h);
        window.clear();

        Moon.setPosition(sf::Vector2f (stateStartParams.positionMoon.x - radiusMoon, stateStartParams.positionMoon.y - radiusMoon));

        for(int i=0; i < stateStartParams.positionEarthPoints.size(); i++) {
            Coordinate vectorDistance = stateStartParams.positionEarthPoints[0] - stateStartParams.positionEarthPoints[i];
            earthCircles[i].setPosition(sf::Vector2f( 20 *vectorDistance.x + (stateStartParams.positionEarthPoints[i].x -  radiusEarth),
                                                     20 *  vectorDistance.y + (stateStartParams.positionEarthPoints[i].y - radiusEarth)));
            window.draw(earthCircles[i]);
        }

        window.draw(earthCircles[0]);
        window.draw(Moon);
        window.display();
        N-=1;

        sf::Event event;
        while (window.pollEvent(event)) {
            switch (event.type) {
                case sf::Event::Closed:
                    window.close();
                    return 0;
            }
        }
        std::cout<<stateStartParams<<std::endl;
        std::cout<<"-------------------------------\n";
    }
    std::cout<<stateStartParams;
    return 0;
}
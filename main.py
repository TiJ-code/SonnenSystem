from vpython import *
from collections import namedtuple
import argparse
import json

# argument parsing
parser = argparse.ArgumentParser(description="Sonnensystem-Simulator")
parser.add_argument("-t", "--time", type=float, default=0, dest="time",
                    help="Wert der zu simulierenden Zeit des Sonnensystems. 0 = Unendlich(Default: 0)")
parser.add_argument("--dt", type=float, default=0.01, dest="dt",
                    help="Zeitschritteeinheit. (Default: 0.01)")
parser.add_argument("--scale", type=float, default=1000, dest="scale",
                    help="Radien skalieren für bessere Sichtbarkeit. Betrifft nur visuelles nicht die Kollisionen. (Default: 1000)")
parser.add_argument("--rate", type=int, default=10000, dest="rate",
                    help="Ausführungen pro Sekunde (Default: 10000)")
parser.add_argument("--configfile", type=str, default="config.json", dest="configfile",
                    help="Pfad zur Config. (Default: config.json)")
parser.add_argument("--useconfig", action="store_true", default=True, dest="useconfig",
                    help="[True] Um die Config zu verwenden. [False] Für nicht. (Default: False)")
parser.add_argument("--integrator", type=str, default="euler", dest="integrator",
                    help="Der zu nutzende Integrator. Optionen: euler, verlet, rk4 (Default: euler)")
parser.add_argument("--endPos", action="store_true", default=False, dest="printEndPos",
                    help="[True] Die End-Position jedes Planeten wird ausgegeben. (Default: False)")
parser.add_argument("--checkEndPos", action="store_true", default=False, dest="checkEndPos",
                    help="[True] Die End-Position jedes Planeten wir d mit der echten End-Position verglichen. (Default: False)")

### containerVector
conVec = namedtuple("conVec", "x y")

### Funktionen
# Gravitationsbedingte Beschleunigung für Euler und Verlet
def gravitational_acc(position):
    sum_acc = vector(0,0,0)
    # Die Geschwindigkeit anhand der anderen Körper berechnen
    for body in bodies:
        # Entfernung zum anderen Körper
        r = mag(position - body.position)
        # Überspringe den Körper wenn er er selbst ist
        if r < body.radius:
            continue
        # Die Länge der Kraft
        acc = G * body.mass / r**2
        # Einheitsvektor der Kraft
        dir = norm(body.position - position)
        # Kraftvektor
        acc = acc * dir
        # Hinzufügen des Vektors zur Summe
        sum_acc += acc
    return sum_acc

# Gravitationsbedingte Beschleunigung für Runge-Kutta
def gravitational_acc_runge(xv):
    sum_acc = vector(0,0,0)
    # Die Geschwindigkeit anhand der anderen Körper berechnen
    for body in bodies:
        # Entfernung zum anderen Körper
        r = mag(xv.x - body.temp_position)
        # Überspringe den Körper wenn er er selbst ist
        if r < body.radius:
            continue
        # Die Länge der Kraft
        acc = G * body.mass / r ** 2
        # Einheitsvektor der Kraft
        dir = norm(body.position - xv.x)
        # Kraftvektor
        acc = acc * dir
        # Hinzufügen des Vektors zur Summe
        sum_acc += acc
    return conVec(xv.y, sum_acc)

### Integratoren
def Euler():
    # Beschleunigung und Geschwindigkeits Berechnungen
    for body in bodies:
        body.acc = gravitational_acc(body.position)
        body.velocity += dt * body.acc

    # Positionen berechnen
    for body in bodies:
        body.position += dt * body.velocity

        # Positionen der Grafiken updaten
        body.sphere.pos = body.position
        body.label.pos = body.position


# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta():
    # berechnen von k1
    for body in bodies:
        body.k = [0]
        body.temp_position = body.position
        body.xv = conVec(body.temp_position, body.velocity)
        temp_k1 = gravitational_acc_runge(body.xv)
        k1 = conVec(temp_k1.x * dt, temp_k1.y * dt)
        body.k.append(k1)

    # verschieben der temp_pos in Richtung k1
    for body in bodies:
        body.temp_position = body.position + body.k[1].x/2

    # berechnen von k2
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k[1].x/2, body.xv.y + body.k[1].y/2)
        temp_k2 = gravitational_acc_runge(temp_xv)
        k2 = conVec(temp_k2.x * dt, temp_k2.y * dt)
        body.k.append(k2)

    # verschieben der temp_pos in Richtung k2
    for body in bodies:
        body.temp_position = body.position + body.k[2].x/2

    # berechnen von k3
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k2[2].x/2, body.xv.y + body.k[2].y/2)
        temp_k3 = gravitational_acc_runge(temp_xv)
        k3 = conVec(temp_k3.x * dt, temp_k3.y * dt)
        body.k.append(k3)

    # verschieben der temp_pos in Richtung k3
    for body in bodies:
        body.temp_position = body.position + body.k[3].x

    # berechnen von k4
    for body in bodies:
        temp_xv = conVec(body.xv.x + body.k[3].x, body.xv.y + body.k[3].y)
        temp_k4 = gravitational_acc_runge(temp_xv)
        k4 = conVec(temp_k4.x * dt, temp_k4 * dt)
        body.k.append(k4)

    # berechnen der Summe und bewegen der Körper
    for body in bodies:
        body.position += 1/6 * (body.k[1].x + 2*body.k[2].x + 2*body.k[3].x + body.k[4].x)
        body.velocity += 1/6 * (body.k[1].y + 2*body.k[2].y + 2*body.k[3].y + body.k[4].y)

        body.sphere.pos = body.position
        body.label.pos = body.position

def Verlet():
    # Positions Berechnung
    for body in bodies:
        body.position += body.velocity * dt + body.acc/2 * dt**2
        body.sphere.pos = body.position
        body.label.pos = body.position
    # Beschleunigung und Geschwindigkeits Berechnungen
    for body in bodies:
        temp_acc = gravitational_acc(body.position)
        body.velocity += dt / 2 * (body.acc + temp_acc)
        body.acc = temp_acc


# Parsing cmd Argumente
args = parser.parse_args()

# Laden der Config
with open(args.configfile, "r") as configfile:
    config = json.load(configfile)

# Konfiguration aus der Config
if args.useconfig:
    dt = config[1]["dt"]
    scale_factor = config[1]["scale_factor"]
    end_time = config[1]["time"]
    integrator = config[1]["integrator"]
else:
    dt = args.dt
    scale_factor = args.scale
    end_time = args.time
    integrator = args.integrator

# Integrator zuweisen
if integrator == "euler":
    integrator = Euler
elif integrator == "rk4":
    integrator = Runge_Kutta
elif integrator == "verlet":
    integrator = Verlet

# EINHEIT:
# Masse: Sonnenmasse
# Länge: Astronomische Einheit
# Zeit: Tage
# G = 4pi^2*AU^3/(M * 365.25) => G = 4*pi^2/365.25^2
# G = 6.674e-11
# scale_factor = 1000
# dt = 0.01

G = 2.9592e-04 # Das 400tsd-fache das der sonst zu gering ist
AU = 1.5e11
M = 2e30
time = 0

# Liste aller Körper der Simulation
bodies = []


class Body():
    def __init__(self, mass=1, radius=1, velocity=vector(0, 0, 0), position=vector(0, 0, 0), color=color.white,
                 trail=True, name="Body", scale=True, index=0):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.temp_position = vector(0, 0, 0)
        self.k = []
        self.xv = conVec(0, 0)
        self.sum_forces = vector(0, 0, 0)
        self.color = color
        self.radius = radius
        # Startwert für die Beschleunigung beim nutzen der Verlet Funktion
        self.acc = gravitational_acc(self.position)
        self.name = name
        self.label = label(pos=self.position, text=self.name, height=10)
        self.index = index
        if scale:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius * scale_factor,
                                 make_trail=trail, retain=200, index=self.index)
        else:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius, make_trail=trail, retain=200, index=self.index)

    # Alternativfunktion zur Euler Integration
    def updateVerlet(self):
        self.forces = []
        self.sum_force = vector(0, 0, 0)
        self.gravitational_force()

        # Kräfte zusammenrechnen
        for force in self.forces:
            self.sum_force += force

        self.position += self.velocity * dt + self.acc / 2 * dt ** 2
        self.velocity += dt / 2 * (self.acc + self.sum_force / self.mass)
        self.acc = self.sum_force / self.mass

        self.sphere.pos = self.position
        self.label.pos = self.position


def color_to_vector(color_list):
    return vector(color_list[0] / 255, color_list[1] / 255, color_list[2] / 255)


for body in config[0]:
    print(body["name"])
    bodies.append(Body(
        name=body["name"],
        mass=body["mass"],
        radius=body["radius"],
        position=vector(body["position"][0], body["position"][1], body["position"][2]),
        velocity=vector(body["velocity"][0], body["velocity"][1], body["velocity"][2]),
        trail=body["trail"],
        color=color_to_vector(body["color"]),
        scale=body["scale"],
        index=len(bodies)
    ))

# erstellt alle Körper der Simulation

time_label = label(pos=vector(75, 350, 0), pixel_pos=True, text="Zeit: " + str(time / 365) + " Jahre")

### Info Box
info_label = label(pos=vector(20, scene.height/3, 0), pixel_pos=True, box=False, align='left', text="")

def onClick(e):
    obj = scene.mouse.pick # Sphere (Planet) auswählen
    if(obj != None):
        # jede Sphere (Planet) hat einen Index, welcher die Position in der Liste aller Körper angibt
        body = bodies[obj.index]
        # Infos über einen Planeten können nicht aus der Klasse kommen
        info_label.text = '<b><i>'+ body.name +'</i></b>\n<i>Masse:</i> '+ str(body.mass) +' Mø\n'
        # TODO: Mehr informationen über die Planeten hinzufügen & Einheiten umrechnen

scene.bind('click', onClick)
###

# Führt die Update Methode jedes Körpers aus
if end_time > 0:
    for epoch in range(int(end_time / dt)):
        rate(args.rate)

        integrator()

        time = epoch * dt
        time_label.text = "Zeit: {:.2f} Jahre".format(time / 365)
    if args.printEndPos:
        for body in bodies:
            print(f"{body.name}: {body.position}")
    if args.checkEndPos:
        error_sum = 0
        for body in bodies:
            end_pos = config[0][body.index]["end_position"]
            error = mag(body.position - end_pos)
            error_sum += error
            print(f"{body.name}: {error} AU")
        print(f"Total error: {error_sum}")
        print(f"dt: {dt}")
        print(f"Integrator: {args.integrator}")

else:
    while True:
        rate(args.rate)

        integrator()

        time_label.text = "Zeit: {:.2f} Jahre".format(time / 365)
        time += dt
from vpython import *
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
parser.add_argument("--useconfig", action="store_true", default=False, dest="useconfig",
                    help="[True] Um die Config zu verwenden. [False] Für nicht. (Default: False)")
parser.add_argument("--integrator", type=str, default="euler", dest="integrator",
                    help="Der zu nutzende Integrator. Optionen: euler, verlet, rk4 (Default: euler)")


### Funktionen
def gravitational_force(bodyself):
    bodyself.sum_forces = vector(0, 0, 0)
    # Berechnen der einwirkenden Gravitation anderer Körper
    for body in bodies:
        # Entfernung zum anderen Körper
        r = mag(bodyself.position - body.position)
        # Überspringe den Körper wenn er er selbst ist
        if r < bodyself.radius + body.radius:
            continue
        # Die Länge der Kraft
        force = G * bodyself.mass * body.mass / r ** 2
        # Einheitsvektor der Kraft
        dir = norm(body.position - bodyself.position)
        # Kraftvektor
        force *= dir
        # Kraftvektor auf die Summenkraft addieren
        bodyself.sum_forces += force


### Integratoren
def Euler():
    # Beschleunigung und Geschwindigkeits Berechnungen
    for body in bodies:
        gravitational_force(body)

        body.acc = body.sum_forces / body.mass
        body.velocity += dt * body.acc

    # Positionen berechnen
    for body in bodies:
        body.position += dt * body.velocity

        # Positionen der Grafiken updaten
        body.sphere.pos = body.position
        body.label.pos = body.position


# Runge-Kutta 4 (RK4) Integrator
def Runge_Kutta():
    pass

def Verlet():
    # Positions Berechnung
    for body in bodies:
        body.position += body.velocity * dt + body.acc / 2 * dt ** 2
        body.sphere.pos = body.position
        body.label.pos = body.position
    # Beschleunigung und Geschwindigkeits Berechnungen
    for body in bodies:
        gravitational_force(body)
        body.velocity += dt / 2 * (body.acc + body.sum_forces / body.mass)
        body.acc = body.sum_forces / body.mass


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
                 trail=True, name="Body", scale=True):
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.sum_forces = vector(0, 0, 0)
        self.color = color
        self.radius = radius
        gravitational_force(self)
        # Startwert für die Beschleunigung beim nutzen der Verlet Funktion
        self.acc = self.sum_forces / self.mass
        self.name = name
        self.label = label(pos=self.position, text=self.name, height=10)
        if scale:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius * scale_factor,
                                 make_trail=trail, retain=5000)
        else:
            self.sphere = sphere(pos=self.position, color=self.color, radius=self.radius, make_trail=trail, retain=5000)

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
        scale=body["scale"]
    ))

# erstellt alle Körper der Simulation

time_label = label(pos=vector(75, 350, 0), pixel_pos=True, text="Zeit: " + str(time / 365) + " Jahre")

# Führt die Update Methode jedes Körpers aus
if end_time > 0:
    for epoch in range(int(end_time / dt)):
        rate(args.rate)

        integrator()

        time = epoch * dt
        time_label.text = "Zeit: {:.2f} Jahre".format(time / 365)


else:
    while True:
        rate(args.rate)

        integrator()

        time_label.text = "Zeit: {:.2f} Jahre".format(time / 365)
        time += dt
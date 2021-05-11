#!/usr/bin/ruby

class BinaryACNode
  attr_reader :matches, :children, :parent
  attr_accessor :index, :failure

  def initialize(parent = nil)
    @index = nil
    @is_root = parent.nil?
    @children = {}
    @parent = if root? then self else parent end
    @failure = if root? then self else nil end
    @matches = []
  end

  def add_child(bit)
    @children[bit] ||= BinaryACNode.new(self)
  end

  def add_match(index, match)
    @matches << [index, match]
  end

  def find_next(bit)
    @children[bit] || if root? then self else failure.find_next bit end
  end

  def root?
    @is_root
  end

  def final?
    !matches.empty?
  end
end

class BinaryACMachine
  attr_reader :root, :current_node

  def initialize(terms, names)
    @root = BinaryACNode.new
    @index2node = []  # Set by method attach_index
    @names = names
    @current_node = @root

    make_goto terms
    make_failure
    attach_index
  end

  def dump_det_wfa(fh)
    fh.puts @index2node.size
    @index2node.each do |n|
      c0 = n.find_next(0)
      c1 = n.find_next(1)
      final = !n.final?
      fh.puts "#{n.index}#{final ? "*" : ""}\t#{c0.index}\t#{c1.index}"
    end
  end

  def reset
    @current_node = @root
  end

  def step(char)
    @current_node = @current_node.find_next(char)
    @current_node.matches
  end

  private

  def make_goto(terms)
    terms.each.with_index do |term, index|
      term.unpack("b*")[0].split("").map(&:to_i).inject(@root) { |node, bit|
        node.add_child bit
      }.add_match(index, term)
    end
  end

  def make_failure
    que = @root.children.to_a
    until que.empty?
      bit, node = que.shift
      node.children.to_a.each do |n|
        que.push n
      end

      node.failure = if node.parent.root?
          @root
        else
          node.parent.failure.find_next bit
        end

      node.failure.matches.each do |m|
        node.add_match *m
      end
    end
  end

  def attach_index
    index = 0
    que = [@root]
    until que.empty?
      node = que.shift
      next unless node.index.nil?
      node.children.each do |_, child|
        que.push child
      end
      node.index = index
      @index2node[node.index] = node
      index += 1
    end
  end
end

src = DATA.lines.map(&:strip).to_a
acm = BinaryACMachine.new src, src
acm.dump_det_wfa($stdout)

"quick fox jumps over the lazy dog".unpack("b*")[0].split("").map(&:to_i).each_with_index do |bit, index|
  final = !acm.step(bit).empty?
  print (final ? 1 : 0) if index % 8 == 7
end
puts

#src = %w!ab bc bab d abcde!
#acm = BinaryACMachine.new src, src
#acm.dump_det_wfa($stdout)

#names = []
#bins = []
#open("min20.db") do |fh|
#  fh.each do |line|
#    next unless line =~ /^([^:]+):([0-9a-f]+)$/
#    names.push $1
#    bins.push [$2].pack("H*")
#  end
#end
#acm = BinaryACMachine.new bins, names
#acm.dump_det_wfa
#$stderr.puts names[-1]


__END__
come
get
give
go
keep
let
make
put
seem
take
be
do
have
say
see
send
may
will
about
across
after
against
among
at
before
between
by
down
from
in
off
on
over
through
to
under
up
with
as
for
of
till
than
a
the
all
any
every
no
other
some
such
that
this
I
he
you
who
and
because
but
or
if
though
while
how
when
where
why
again
ever
far
forward
here
near
now
out
still
then
there
together
well
almost
enough
even
little
much
not
only
quite
so
very
tomorrow
yesterday
north
south
east
west
please
yes
account
act
addition
adjustment
advertisement
agreement
air
amount
amusement
animal
answer
apparatus
approval
argument
art
attack
attempt
attention
attraction
authority
back
balance
base
behaviour
belief
birth
bit
bite
blood
blow
body
brass
bread
breath
brother
building
burn
burst
business
butter
canvas
care
cause
chalk
chance
change
cloth
coal
colour
comfort
committee
company
comparison
competition
condition
connection
control
cook
copper
copy
cork
cotton
cough
country
cover
crack
credit
crime
crush
cry
current
curve
damage
danger
daughter
day
death
debt
decision
degree
design
desire
destruction
detail
development
digestion
direction
discovery
discussion
disease
disgust
distance
distribution
division
doubt
drink
driving
dust
earth
edge
education
effect
end
error
event
example
exchange
existence
expansion
experience
expert
fact
fall
family
father
fear
feeling
fiction
field
fight
fire
flame
flight
flower
fold
food
force
form
friend
front
fruit
glass
gold
government
grain
grass
grip
group
growth
guide
harbour
harmony
hate
hearing
heat
help
history
hole
hope
hour
humour
ice
idea
impulse
increase
industry
ink
insect
instrument
insurance
interest
invention
iron
jelly
join
journey
judge
jump
kick
kiss
knowledge
land
language
laugh
law
lead
learning
leather
letter
level
lift
light
limit
linen
liquid
list
look
loss
love
machine
man
manager
mark
market
mass
meal
measure
meat
meeting
memory
metal
middle
milk
mind
mine
minute
mist
money
month
morning
mother
motion
mountain
move
music
name
nation
need
news
night
noise
note
number
observation
offer
oil
operation
opinion
order
organization
ornament
owner
page
pain
paint
paper
part
paste
payment
peace
person
place
plant
play
pleasure
point
poison
polish
porter
position
powder
power
price
print
process
produce
profit
property
prose
protest
pull
punishment
purpose
push
quality
question
rain
range
rate
ray
reaction
reading
reason
record
regret
relation
religion
representative
request
respect
rest
reward
rhythm
rice
river
road
roll
room
rub
rule
run
salt
sand
scale
science
sea
seat
secretary
selection
self
sense
servant
sex
shade
shake
shame
shock
side
sign
silk
silver
sister
size
sky
sleep
slip
slope
smash
smell
smile
smoke
sneeze
snow
soap
society
son
song
sort
sound
soup
space
stage
start
statement
steam
steel
step
stitch
stone
stop
story
stretch
structure
substance
sugar
suggestion
summer
support
surprise
swim
system
talk
taste
tax
teaching
tendency
test
theory
thing
thought
thunder
time
tin
top
touch
trade
transport
trick
trouble
turn
twist
unit
use
value
verse
vessel
view
voice
walk
war
wash
waste
water
wave
wax
way
weather
week
weight
wind
wine
winter
woman
wood
wool
word
work
wound
writing
year
angle
ant
apple
arch
arm
army
baby
bag
ball
band
basin
basket
bath
bed
bee
bell
berry
bird
blade
board
boat
bone
book
boot
bottle
box
boy
brain
brake
branch
brick
bridge
brush
bucket
bulb
button
cake
camera
card
cart
carriage
cat
chain
cheese
chest
chin
church
circle
clock
cloud
coat
collar
comb
cord
cow
cup
curtain
cushion
dog
door
drain
drawer
dress
drop
ear
egg
engine
eye
face
farm
feather
finger
fish
flag
floor
fly
foot
fork
fowl
frame
garden
girl
glove
goat
gun
hair
hammer
hand
hat
head
heart
hook
horn
horse
hospital
house
island
jewel
kettle
key
knee
knife
knot
leaf
leg
library
line
lip
lock
map
match
monkey
moon
mouth
muscle
nail
neck
needle
nerve
net
nose
nut
office
orange
oven
parcel
pen
pencil
picture
pig
pin
pipe
plane
plate
plough
pocket
pot
potato
prison
pump
rail
rat
receipt
ring
rod
roof
root
sail
school
scissors
screw
seed
sheep
shelf
ship
shirt
shoe
skin
skirt
snake
sock
spade
sponge
spoon
spring
square
stamp
star
station
stem
stick
stocking
stomach
store
street
sun
table
tail
thread
throat
thumb
ticket
toe
tongue
tooth
town
train
tray
tree
trousers
umbrella
wall
watch
wheel
whip
whistle
window
wing
wire
worm
able
acid
angry
automatic
beautiful
black
boiling
bright
broken
brown
cheap
chemical
chief
clean
clear
common
complex
conscious
cut
deep
dependent
early
elastic
electric
equal
fat
fertile
first
fixed
flat
free
frequent
full
general
good
great
grey
hanging
happy
hard
healthy
high
hollow
important
kind
like
living
long
male
married
material
medical
military
natural
necessary
new
normal
open
parallel
past
physical
political
poor
possible
present
private
probable
quick
quiet
ready
red
regular
responsible
right
round
same
second
separate
serious
sharp
smooth
sticky
stiff
straight
strong
sudden
sweet
tall
thick
tight
tired
true
violent
waiting
warm
wet
wide
wise
yellow
young
awake
bad
bent
bitter
blue
certain
cold
complete
cruel
dark
dead
dear
delicate
different
dirty
dry
false
feeble
female
foolish
future
green
ill
last
late
left
loose
loud
low
mixed
narrow
old
opposite
public
rough
sad
safe
secret
short
shut
simple
slow
small
soft
solid
special
strange
thin
white
wrong

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
end

class BinaryACMachine
  attr_reader :root

  def initialize(terms, names)
    @root = BinaryACNode.new
    @index2node = []  # Set by method attach_index
    @names = names

    make_goto terms
    make_failure
    attach_index
  end

  def dump_det_wfa
    table = @index2node.map do |n|
      c0 = n.find_next(0)
      c1 = n.find_next(1)
      g0 = if c0.matches.empty? then 0 else 1 end
      g1 = if c1.matches.empty? then 0 else 1 end
      [n.index, c0.index, c1.index, g0, g1, n.matches]
    end
    table.each do |index, child0, child1, g0, g1, matches|
      puts "{#{index}, #{child0}, #{child1}, #{g0}, #{g1}},\t// #{matches.map { |m| @names[m[0]] }.join(", ")}"
    end
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

src = %w!ab bc bab d abcde!
acm = BinaryACMachine.new src, src
acm.dump_det_wfa

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

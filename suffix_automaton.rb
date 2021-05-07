#!/usr/bin/ruby

require "stringio"

class State
  attr_accessor :len, :link, :next

  def initialize(len: 0, link: 0)
    @len = len
    @link = link
    @next = {}
  end
end

# Thanks to: https://w.atwiki.jp/uwicoder/pages/2842.html
class SuffixAutomaton
  def initialize
    @last = 0
    @st = [
      State.new(len: 0, link: -1),
    ]
  end

  def extend(c)
    # FIXME: Make it more Ruby-ish
    cur = @st.size
    @st.push State.new
    @st[cur].len = @st[@last].len + 1
    p = @last
    while p != -1 and not @st[p].next.key?(c)
      @st[p].next[c] = cur
      p = @st[p].link
    end
    if p == -1
      @st[cur].link = 0
    else
      q = @st[p].next[c]
      if @st[p].len + 1 == @st[q].len
        @st[cur].link = q
      else
        clone = @st.size
        @st.push State.new
        @st[clone].len = @st[p].len + 1
        @st[clone].next = @st[q].next.dup
        @st[clone].link = @st[q].link
        while p != -1 && @st[p].next[c] == q
          @st[p].next[c] = clone
          p = @st[p].link
        end
        @st[q].link = @st[cur].link = clone
      end
    end
    @last = cur
  end

  def check(input)
    cur = 0
    input.each do |i|
      n = @st[cur].next[i]
      if n.nil?
        break false
      else
        cur = n
        true
      end
    end
  end

  def dump_det_wfa(fh)
    fail_st = @st.size
    table = @st.map.with_index do |s, si|
      c0 = s.next[0] || fail_st
      c1 = s.next[1] || fail_st
      [si, c0, c1, 1]
    end
    table.push [fail_st, fail_st, fail_st, 0]

    table.each do |index, child0, child1, marker|
      fh.puts "{#{index}, #{child0}, #{child1}, #{marker}},"
    end
  end

  def to_dot(print_next: true, print_link: false)
    sio = StringIO.new
    sio.puts "digraph{"
    sio.puts "graph[rankdir=LR];"
    sio.puts "node[shape=circle];"
    @st.each_with_index do |s, si|
      if print_link and s.link != -1
        sio.puts "\"#{si}\"->\"#{s.link}\"[style=dashed];"
      end
      if print_next
        s.next.each do |nc, ns|
          sio.puts "\"#{si}\"->\"#{ns}\"[label=\"#{nc}\"];"
        end
      end
    end
    sio.puts "}"
    sio.rewind
    sio.read
  end
end

raise "Usage: $0 OUT-FILE-PREFIX" unless ARGV.size == 1
outfile_prefix = ARGV[0]

Random.srand(0)
src = Random.bytes(1000).unpack("b*")[0].split("").map(&:to_i)

m = SuffixAutomaton.new
src.each do |b|
  m.extend b
end
#puts m.to_dot

open("#{outfile_prefix}_automaton.inc", "w") do |fh|
  m.dump_det_wfa fh
end

size = 50 * 8
from = Random.rand(0..src.size - size)
to = from + size
raise "hoge" unless m.check(src[from...to])
open("#{outfile_prefix}_input_valid.inc", "w") do |fh|
  fh.puts src[from...to].map { |i| i != 0 }.join(", ")
end

while true
  cand = Random.bytes(size / 8).unpack("b*")[0].split("").map(&:to_i)
  next if m.check cand
  open("#{outfile_prefix}_input_invalid.inc", "w") do |fh|
    fh.puts cand.map { |i| i != 0 }.join(", ")
  end
  break
end

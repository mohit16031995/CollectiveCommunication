#include <assert.h>
#include <fstream>
#include <graphviz/cgraph.h>
#include <inttypes.h>
#include <iostream>
#include <map>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../../cmdline.h"

#define noorder 1
#define DEBUG 0

const uint8_t NODE_TYPE_HOST = 1;
const uint8_t NODE_TYPE_SWITCH = 2;

class TGNode {

  friend class TopoGraph; // no need to write get/set functions for everything
                          // (c++ rocks, try this in java...)

private:
#ifdef DEBUG
  std::string name;
#endif
  std::vector<TGNode *>
      ports; // ports[outport] = ptr to node connected to this port
  std::vector<uint32_t> linkids; // linkids[port] = the linkid of the link that
                                 // is connected to port
  std::map</*target*/ uint32_t, uint32_t /*port*/> routinginfo;
  uint8_t type;
  uint32_t id;

  size_t node_name_len;
  char *node_name;
};

class Mapping {};

class TopoGraph {

private:
  std::vector<TGNode *> nodes; // nodes[id] = ptr to node with nodeid = id
  std::map<std::string, uint32_t> nodename2nodeid;
  std::map</*linkid*/ uint32_t,
           /*(nodeid, port)*/ std::pair<uint32_t, uint32_t>>
      linkid2node_port;
  std::map<uint32_t, uint32_t>
      rank2nodeid; // switches also have a nodeid, so we need this map
  std::map<uint32_t, uint32_t> nodeid2rank;

  std::map<uint32_t, std::string> rank2nodename;
  std::map<std::string, uint32_t> nodename2rank;
  uint32_t last_free_rank;

public:
  TopoGraph(char *filename, char *mapping_filename, bool print)
      : last_free_rank(0) {
    /* This needs to be before the parsing of dotfile */
    if (mapping_filename != NULL)
      this->parse_mapping_file(mapping_filename);

    this->parse_dotfile_agraph(filename);
  }

  ~TopoGraph() {
    for (auto node : nodes) {
      delete node;
    }
  }

  inline uint32_t get_free_rank() {

    std::map<uint32_t, std::string>::iterator it;

    do {
      it = rank2nodename.find(last_free_rank);
      last_free_rank++;
    } while (it != rank2nodename.end());

    return last_free_rank - 1;
  }

  inline uint32_t create_node(std::string nodename, uint8_t nodetype) {

    /*
        Creates a new node. Nodename must be a unique name, returns the
        id of the created node.
    */

    TGNode *n = new TGNode;

    n->type = nodetype;
    n->id = this->nodes.size();
    n->node_name = NULL;
    n->node_name_len = 0;

#ifdef DEBUG
    n->name = nodename;
#endif
    this->nodes.push_back(n);
    this->nodename2nodeid.insert(std::make_pair(nodename, n->id));

    return n->id;
  }

  int get_size() { return nodes.size(); }

  int get_outports(int i) { return nodes[i]->ports.size(); }

  char *get_node_name(int id, uint32_t *len) {
    *len = nodes[id]->node_name_len;
    return nodes[id]->node_name;
  }
  char *get_node_name(int id) { return nodes[id]->node_name; }

  std::string get_nodename(int id) {
#ifdef DEBUG
    return nodes[id]->name;
#else
    return "debug disabled";
#endif
  }
  inline uint32_t create_link(std::string from, std::string to) {

    /*
        Creates a new directed link between the nodes "from" and "to".
        Returns the link id.
    */

    uint32_t f, t;

    f = (this->nodename2nodeid.find(from))->second;
    t = (this->nodename2nodeid.find(to))->second;

    (this->nodes[f])->linkids.push_back(this->linkid2node_port.size());
    this->linkid2node_port.insert(std::make_pair(
        this->linkid2node_port.size(),
        std::make_pair((this->nodes[f])->id, (this->nodes[f])->ports.size())));
    (this->nodes[f])->ports.push_back(this->nodes[t]);

    return (this->linkid2node_port.size() - 1);
  }

  inline void add_routing_info(uint32_t dest, uint32_t link) {

    /*
        Adds the target with node-id "dest" to the link with the id
        "link". This means if there is a packet for dest at the node
        at tail(link), this packet will use link to reach the next hop.
    */

    std::map</*linkid*/ uint32_t,
             /*(nodeid, port)*/ std::pair<uint32_t, uint32_t>>::iterator
        linkid2node_port_iter;
    linkid2node_port_iter = linkid2node_port.find(link);
    assert(linkid2node_port_iter != linkid2node_port.end());
    uint32_t node = linkid2node_port_iter->second.first;
    uint32_t port = linkid2node_port_iter->second.second;
    nodes[node]->routinginfo.insert(std::make_pair(dest, port));
  }
  /*
  inline void get_nodename(Agnode_t * node){
      Agedge_t* link = agfstout(graph, node_from);

  }*/

  void parse_dotfile_agraph(char *filename) {

    FILE *fd = fopen(filename, "r");
    if (fd == NULL) {
      fprintf(stderr, "Couldn't open %s for reading!\n", filename);
      exit(EXIT_FAILURE);
    }

    Agraph_t *graph = agread(fd, NULL);
    if (graph == NULL) {
      fprintf(stderr, "Couldn't parse graph data!\n");
      fclose(fd);
      exit(EXIT_FAILURE);
    }
    fclose(fd);

    // iterate over graphs node's and add them
    Agnode_t *node = agfstnode(graph);

    while (node != NULL) {
      std::string nodename;
      uint8_t nodetype;
      uint32_t rank;
      nodename = (std::string)agnameof(node);
      if (nodename.find('H', 0) == 0)
        nodetype = NODE_TYPE_HOST;
      else
        nodetype = NODE_TYPE_SWITCH;
      uint32_t nodeid = create_node(nodename, nodetype);
      if (nodetype == NODE_TYPE_HOST) {
        // printf("[NET] adding host: %s; nodeid: %i; rank: %lu\n",
        // nodename.c_str(), nodeid, rank2nodeid.size());

        // try to get rank->nodeid mapping from provided mapping file
        Agedge_t *link = agfstout(graph, node);
        char *hostname = agget(link, ((char *)"hostname"));
        if (hostname != NULL) {
          std::map<std::string, uint32_t>::iterator it =
              nodename2rank.find(hostname);

          if (it == nodename2rank.end())
            rank = get_free_rank();
          else
            rank = it->second;
        } else {
          rank = get_free_rank();
        }
        // printf("[NET] rank %i mapped to %s\n", rank, hostname);
        rank2nodeid.insert(std::make_pair(rank, nodeid));
      } else {
        // printf("[NET] adding switch: %s; nodeid: %i;\n", nodename.c_str(),
        // nodeid);
      }
      node = agnxtnode(graph, node);
    }
    // fprintf(stderr, "rank2nodeid.size(): %d\n", (int)rank2nodeid.size());
    // now we added all the nodes, so we can add the links between them now
    Agnode_t *node_from = agfstnode(graph);

    while (node_from != NULL) {

      std::string node_from_name;

      node_from_name = (std::string)agnameof(node_from);
      Agedge_t *link = agfstout(graph, node_from);

      while (link != NULL) {
        Agnode_t *node_to = aghead(link);
        std::string node_to_name = (std::string)agnameof(node_to);
        uint32_t newlinkid = create_link(node_from_name, node_to_name);
        // add the routinginfo for this link
        char *comment = agget(link, ((char *)"comment"));
        // if the comment is * we leave the routinginfo empty - that means there
        // is only one connection (port 0) and everything goes there

        if (strcmp(comment, "*") != 0) {
          char *buffer = (char *)malloc(strlen(comment) * sizeof(char) + 1);
          strcpy(buffer, comment);
          char *result = strtok(buffer, ", \t\n");
          while (result != NULL) {
            uint32_t dest =
                (this->nodename2nodeid.find(std::string(result)))->second;
            add_routing_info(dest, newlinkid);
            result = strtok(NULL, ", \t\n");
          }
          free(buffer);
        }
        link = agnxtout(graph, link);
      }
      node_from = agnxtnode(graph, node_from);
    }
    agclose(graph);
  }

  uint32_t get_id(uint32_t rank) {
    std::map<uint32_t, uint32_t>::iterator r2nid_it;
    r2nid_it = rank2nodeid.find(rank);
    assert(r2nid_it != rank2nodeid.end());
    return r2nid_it->second;
  }

  uint32_t next_port(uint32_t node_id, uint32_t dest_id) {
    uint32_t pos = node_id;

    if (nodes[pos]->ports.size() == 1) {
      // there is only one outgoing link, use it
      return 0;
    } else {
      // there are multiple outgoing ports - find the right one
      std::map</*target*/ uint32_t, uint32_t /*port*/>::iterator nexthop;
      nexthop = nodes[pos]->routinginfo.find(dest_id);
      assert(nexthop != nodes[pos]->routinginfo.end());
      return nexthop->second;
    }
  }

  uint32_t next_hop(uint32_t node_id, uint32_t dest_id) {
    uint32_t np = next_port(node_id, dest_id);
    return nodes[node_id]->ports[np]->id;
  }

  std::vector</*linkids*/ uint32_t> find_route(uint32_t src_rank,
                                               uint32_t dest_rank) {

    std::map<uint32_t, uint32_t>::iterator r2nid_it;
    r2nid_it = rank2nodeid.find(src_rank);
    assert(r2nid_it != rank2nodeid.end());
    uint32_t src_id = r2nid_it->second;
    r2nid_it = rank2nodeid.find(dest_rank);
    assert(r2nid_it != rank2nodeid.end());
    uint32_t dest_id = r2nid_it->second;

    assert(nodes.size() > dest_id);
    std::vector<uint32_t> path;
    uint32_t pos = src_id;
    while (pos != dest_id) {
      assert(nodes.size() > pos);
      assert(nodes[pos]->ports.size() > 0);
      if (nodes[pos]->ports.size() == 1) {
        // there is only one outgoing link, use it
        path.push_back(nodes[pos]->linkids[0]);
        pos = nodes[pos]->ports[0]->id;
      } else {
        // there are multiple outgoing ports - find the right one
        std::map</*target*/ uint32_t, uint32_t /*port*/>::iterator nexthop;
        nexthop = nodes[pos]->routinginfo.find(dest_id);
        assert(nexthop != nodes[pos]->routinginfo.end());
        path.push_back(nodes[pos]->linkids[nexthop->second]);
        pos = nodes[pos]->ports[nexthop->second]->id;
      }
    }
    return path;
  }

  // rank name
  void parse_mapping_file(std::string filename) {
    std::fstream f;
    f.open(filename, std::fstream::in);
    if (!f.is_open())
      return;

    std::string name;
    uint32_t rank;

    while (f >> rank >> name) {
      nodename2rank.insert(std::make_pair(name, rank));
      // printf("rank: %i; name: %s\n", rank, name.c_str());
      rank2nodename.insert(std::make_pair(rank, name));
    }
    f.close();
  }
};

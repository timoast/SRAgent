# import 
import os
import operator
from enum import Enum
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Sequence, TypedDict
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from pydantic import BaseModel
from langchain_openai import ChatOpenAI
from langgraph.graph import START, END, StateGraph
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from Bio import Entrez
## package
from SRAgent.agents.entrez import create_entrez_agent

# state
class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    messages: Annotated[Sequence[BaseMessage], operator.add, "Messages"]
    SRP: Annotated[List[str], operator.add, "SRP accessions"]
    SRX: Annotated[List[str], operator.add, "SRX accessions"]
    SRR: Annotated[List[str], operator.add, "SRR accessions"]
    route: Annotated[str, "Route choice"]
    rounds: Annotated[int, "Rounds of trying to convert accessions"]


# functions
## entrez agent
def invoke_entrez_agent_node(state: GraphState):
    entrez_agent = create_entrez_agent()
    return entrez_agent.invoke({"messages" : state["messages"]})

## accessions extraction
class Acessions(BaseModel):
    srp: List[str]
    srx: List[str]
    srr: List[str]

def create_get_accessions_node():
    model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)
    def invoke_get_accessions_node(state: GraphState):
        """
        Structured data extraction
        """
        message = state["messages"][-1].content
        prompt = "\n".join([
            "Sequence Read Archive accessions:",
            " - SRP: project accession",
            " - SRX: experiment accession",
            " - SRR: run accession",
            "The most important accession are SRX.",
            f"Extract SRA accessions from the following:",
            message
        ])
    
        response = model.with_structured_output(Acessions, strict=True).invoke(prompt)
        return {
            "SRP" : response.srp,
            "SRX" : response.srx,
            "SRR" : response.srr
        }
    return invoke_get_accessions_node

## router
class Choices(Enum):
    CONTINUE = "Continue"
    STOP = "Stop"

class Choice(BaseModel):
    Choice: Choices

def create_router_node():
    """
    Router for the graph
    """
    model = ChatOpenAI(model="gpt-4o-mini", temperature=0)

    def invoke_router(
        state: GraphState
    ) -> Annotated[dict, "Response from the router"]:
        """
        Route the conversation to the appropriate tool based on the current state of the conversation.
        """
        # create prompt
        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", "You determine whether Sequence Read Archive SRX accessions have been obtained from the Entrez ID."),
            ("system", "Here are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            # Add the final question/instruction
            ("human", "Based on the messages above, select STOP if the task is complete or CONTINUE if more information is needed."),
        ])
        formatted_prompt = prompt.format_messages(history=state["messages"][-4:])
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        return {"route": response.Choice.value, "rounds": 1}
    
    return invoke_router

def route_interpret(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    if state["rounds"] >= 2:
        return END
    return "entrez_agent_node" if state["route"] == "Continue" else END

def create_convert_graph():
    """
    Create a graph that converts Entrez IDs & non-SRA accessions to SRA accessions
    """
    workflow = StateGraph(GraphState)
    
    # nodes
    workflow.add_node("entrez_agent_node", invoke_entrez_agent_node)
    workflow.add_node("get_accessions_node", create_get_accessions_node())
    workflow.add_node("router_node", create_router_node())

    # edges
    workflow.add_edge(START, "entrez_agent_node")
    workflow.add_edge("entrez_agent_node", "get_accessions_node")
    workflow.add_edge("get_accessions_node", "router_node")
    workflow.add_conditional_edges("router_node", route_interpret)

    # compile the graph
    graph = workflow.compile()
    return graph

def invoke_convert_graph(
    state: GraphState, 
    graph: StateGraph,
    to_return: List[str] = ["SRX"]
) -> Annotated[dict, "Response from the graph"]:
    """
    Invoke the graph to convert Entrez IDs & non-SRA accessions to SRA accessions
    """
    response = graph.invoke(state)
    filtered_response = {key: response[key] for key in to_return}
    return filtered_response

# main
if __name__ == "__main__":
    from functools import partial

    #-- setup --#
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    #-- graph --#
    entrez_id = "34748561"
    msg = f"Obtain all SRX accessions for the Entrez ID {entrez_id}"
    input = {"messages" : [HumanMessage(content=msg)]}
    graph = create_convert_graph()
    #for step in graph.stream(input, config={"max_concurrency" : 3, "recursion_limit": 30}):
    #    print(step)

    #invoke_convert_graph = partial(invoke_convert_graph, graph=graph)
    #print(invoke_convert_graph(input))
    

    #-- nodes --#
    state = {
        "entrez_id" : "34748561",
        "messages" : []
    }
    # entrez_agent_node = create_entrez_agent_node(state)

    state = {
        "messages" : [
            HumanMessage(
                content=" - **Study Accession**: SRP526682\n- **Sample Accession**: SRS22358714- **Run Accession**: SRR30256648\n**Experiment Accession**: SRX4967528"
            )
        ]
    }
    # create_get_accessions_node(state)

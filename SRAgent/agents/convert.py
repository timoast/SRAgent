# import 
import os
import operator
from enum import Enum
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Sequence, TypedDict, Callable
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from pydantic import BaseModel
from langchain_openai import ChatOpenAI
from langgraph.graph import START, END, StateGraph
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
## package
from SRAgent.agents.entrez import create_entrez_agent

# state
class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    messages: Annotated[Sequence[BaseMessage], operator.add]
    entrez_id: Annotated[str, "Entrez ID"]
    SRP: Annotated[List[str], operator.add]
    SRX: Annotated[List[str], operator.add]
    SRR: Annotated[List[str], operator.add]
    route: Annotated[str, "Route choice"]
    rounds: Annotated[int, operator.add]


# functions
## entrez agent
def invoke_entrez_agent_node(state: GraphState) -> Dict[str, List[str]]:
    entrez_agent = create_entrez_agent()
    response = entrez_agent.invoke({"messages" : state["messages"]})
    return {"messages" : [response["messages"][-1]]}

## accessions extraction
class Acessions(BaseModel):
    srp: List[str]
    srx: List[str]
    srr: List[str]

def create_get_accessions_node() -> Callable:
    model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)
    def invoke_get_accessions_node(state: GraphState):
        """
        Structured data extraction
        """
        # create prompt
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
        # invoke model with structured output
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
    Message: str

def create_router_node() -> Callable:
    """
    Router for the graph
    """
    model = ChatOpenAI(model="gpt-4o", temperature=0)

    def invoke_router(
        state: GraphState
    ) -> Annotated[dict, "Response from the router"]:
        """
        Route the conversation to the appropriate tool based on the current state of the conversation.
        """
        def format_accessions(accessions):
            if not accessions:
                return "No accessions found"
            return ", ".join(accessions)

        accesions = "\n".join([
            " - SRP: " + format_accessions(state["SRP"]),
            " - SRX: " + format_accessions(state["SRX"]),
            " - SRR: " + format_accessions(state["SRR"]),
            "\n"
        ])

        # create prompt
        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", 
                "You determine whether Sequence Read Archive SRX accessions (e.g., SRX123456) have been obtained from the Entrez ID."
                " There should be at least one SRX accession. SRP and SRR accessions are optional."
                " If the accessions have been obtained, select STOP. If more information is needed, select CONTINUE."
                " If more information is needed (CONTINUE), provide one or two sentences of feedback on how to obtain the data (e.g., use esearch instead of efetch)."),
            ("system", "Here are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            ("system", "\nHere are the extracted SRA accessions:\n" + accesions)
        ])
        formatted_prompt = prompt.format_messages(history=state["messages"][-4:])
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        # format the response
        return {"route": response.Choice.value,  "messages": [AIMessage(content=response.Message)], "rounds": 1}
    
    return invoke_router

def route_interpret(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    if state["rounds"] >= 2:
        return END
    return "entrez_agent_node" if state["route"] == "Continue" else END

def create_convert_graph() -> StateGraph:
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
    # filter to just the keys we want to return
    return {key: [response[key]] for key in to_return}

# main
if __name__ == "__main__":
    from functools import partial
    from Bio import Entrez

    #-- setup --#
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    #-- graph --#
    entrez_id = "34748561"
    msg = f"Obtain all SRX accessions for the Entrez ID {entrez_id}"
    input = {"messages" : [HumanMessage(content=msg)]}
    graph = create_convert_graph()
    for step in graph.stream(input, config={"max_concurrency" : 3, "recursion_limit": 30}):
        print(step)

    ## invoke with graph object directly provided
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
